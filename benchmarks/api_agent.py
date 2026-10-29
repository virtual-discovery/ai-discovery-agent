"""Shared Harbor shell agent for API models on the benchmark.

Provider credentials stay on the host. Full request/response and command output
are recorded; large command output is explicitly paginated, never silently lost.
API limits and incomplete responses fail the smoke instead of being scored as a
normal completed episode. This is a benchmark harness, not AI Theorist.
"""
import asyncio
import json
import os
from pathlib import Path
import time
import urllib.error
import urllib.request

from harbor.agents.base import BaseAgent
from benchmarks import durable_shell

ROOT = Path(__file__).resolve().parent.parent
SYSTEM = """You are a scientific computing agent completing an isolated benchmark.
Work in /app using the shell tool. Read the task's files and actually run the
calculations; do not invent tool results. Follow the task's budget and submission
contract. Python, numpy, scipy and Quantum ESPRESSO are installed. No internet is
available inside the task. Finish only after saving the requested artifacts.
Shell commands run durably. A long command returns state=running and a command_id;
use poll_command with that ID until it completes. Do not restart a running
command or the calculation it is polling. poll_command waits up to 240 seconds;
do useful independent work between polls when possible. Complete stdout and
stderr are retained in files and returned after completion, even for commands
that run beyond 120 seconds. Tool output is paginated at 16000 characters. When
has_more is true, use read_output to retrieve subsequent pages; full output is
retained. This paging is not a signal that the command failed. Do not access
private lab or verifier files; use the documented scientific interface.
"""
TOOLS = [
    {"name": "shell", "description": "Execute a shell command in /app in the task container.",
     "parameters": {"type": "object", "properties": {"command": {"type": "string"}},
                    "required": ["command"], "additionalProperties": False}},
    {"name": "poll_command", "description": "Wait up to 240 seconds for an existing shell command; return running status or its complete paginated output. Never restarts it.",
     "parameters": {"type": "object", "properties": {"command_id": {"type": "integer"}},
                    "required": ["command_id"], "additionalProperties": False}},
    {"name": "read_output", "description": "Read the next page of a previous command's complete output.",
     "parameters": {"type": "object", "properties": {"output_id": {"type": "integer"},
                     "offset": {"type": "integer"}}, "required": ["output_id", "offset"],
                    "additionalProperties": False}},
]


def credentials(name):
    # Host-only credentials; never copy them into scientific task containers.
    if os.environ.get(name, '').strip():
        return os.environ[name].strip()
    path = ROOT / '.env'
    for line in path.read_text().splitlines() if path.exists() else []:
        key, sep, value = line.strip().removeprefix('export ').partition('=')
        if sep and key.strip() == name:
            value = value.strip().strip('"').strip("'")
            if value and not value.startswith('YOUR_'):
                return value
    raise RuntimeError(f'Set {name} on the host or in the repository .env')


# Long enough that slow-host SCFs do not turn into dozens of polling calls,
# short enough to keep the 5-minute provider prompt caches warm.
POLL_WAIT_SECONDS = 240
TRANSIENT_HTTP = {408, 429, 500, 502, 503, 504, 520, 521, 522, 523, 524, 529}  # 52x: Cloudflare edge errors
# USD per million tokens. Claude: input, cache read, cache write, output.
# OpenAI: input, cached input, output, then the long-context threshold and
# long-context rates (None when the model has no long-context tier).
CLAUDE_RATES = {'claude-opus-4-6': (5, .5, 6.25, 25), 'claude-sonnet-5': (2, .2, 2.5, 10), 'claude-opus-5-5': (4, .2, 5, 20)}
OPENAI_RATES = {'gpt-5.3-codex': (1.75, .175, 14, None, None),
                'gpt-5.6-terra': (2, .2, 12, 272000, (4, .4, 18)),
                'gpt-6-astra': (10, 1, 50, None, None)}   # confirm the API model id before use
RETRY_DELAYS = (5, 15, 45, 120)


def post(url, headers, payload, retries=None, timeout=240):
    # Bounded retries on transient provider/network errors only; unbilled
    # failed requests are logged in `retries` and never change the transcript.
    data = json.dumps(payload).encode()
    for attempt in range(len(RETRY_DELAYS) + 1):
        request = urllib.request.Request(url, data=data,
                                         headers={'Content-Type': 'application/json', **headers})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except urllib.error.HTTPError as exc:
            # Avoid logging headers, keys or provider echoes in exceptions.
            error = f'Provider HTTP {exc.code}'
            transient = exc.code in TRANSIENT_HTTP
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            error = f'Provider network error {type(exc).__name__}'
            transient = True
        if not transient or attempt == len(RETRY_DELAYS):
            raise RuntimeError(error) from None
        if retries is not None:
            retries.append({'error': error, 'attempt': attempt + 1, 'epoch': time.time()})
        time.sleep(RETRY_DELAYS[attempt])


class PaperBaselineAgent(BaseAgent):
    @staticmethod
    def name():
        return 'paper-baseline'

    def version(self):
        return '0.8.0'

    def __init__(self, *args, max_calls=30, max_cost_usd=5, prompt_cache=True, agent_qe=True,
                 max_output_tokens=12288, truncation_retries=0, **kwargs):
        super().__init__(*args, **kwargs)
        self.max_calls = int(max_calls)
        self.max_cost = float(max_cost_usd)
        self.max_output = int(max_output_tokens)
        # Default: a truncated response fails the episode. The runner allows re-issuing the identical request
        # a bounded number of times (each attempt is billed, recorded and counts as a call).
        self.truncation_retries = int(truncation_retries)
        # Benchmark agent images have no Quantum ESPRESSO; calculations go through the lab service.
        self.system = SYSTEM if str(agent_qe).lower() not in ('false', '0', 'no') else SYSTEM.replace(
            'Python, numpy, scipy and Quantum ESPRESSO are installed.', "Python, numpy and scipy are installed; calculations run through the task's lab service.")
        self.prompt_cache = str(prompt_cache).lower() not in ('false', '0', 'no')
        self.outputs = []
        self.commands = {}

    async def setup(self, environment):
        result = await environment.exec('python3 -c "import numpy, scipy; print(\'runtime-ready\')"', timeout_sec=20)
        if result.return_code:
            raise RuntimeError('Task Python runtime unavailable')

    def save(self, name, value):
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        (self.logs_dir / name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')

    def stop_at_resource_limit(self, context, reason):
        # A bounded agent failing to finish is a benchmark outcome, not a broken
        # evaluator. Let Harbor grade its actual artifacts (usually reward zero).
        # Keep completed=False; never count this as a completed API smoke.
        context.metadata['termination_reason'] = reason
        self.save('usage.json', context.model_dump())

    def page(self, output_id, offset):
        if not isinstance(output_id, int) or not 0 <= output_id < len(self.outputs):
            return {'error': 'Unknown output_id'}
        output = self.outputs[output_id]
        if output is None:
            return {'state': 'running', 'command_id': output_id,
                    'instruction': 'Use poll_command to await completion; full output is retained.'}
        if not isinstance(offset, int) or not 0 <= offset <= len(output):
            return {'error': 'Offset out of bounds'}
        end = min(offset + 16000, len(output))
        return {'output_id': output_id, 'offset': offset, 'next_offset': end,
                'total_characters': len(output), 'has_more': end < len(output),
                'output': output[offset:end]}

    async def execute(self, name, args, environment):
        if name == 'read_output':
            return self.page(args.get('output_id'), args.get('offset'))
        if name == 'poll_command':
            command_id = args.get('command_id')
            if type(command_id) is not int or command_id not in self.commands:
                return {'error': 'Unknown command_id'}
            return await self.poll_command(command_id, environment, POLL_WAIT_SECONDS)
        if name != 'shell' or not isinstance(args.get('command'), str):
            return {'error': 'Invalid tool invocation'}
        output_id = len(self.outputs)
        self.outputs.append(None)
        self.commands[output_id] = await durable_shell.launch(environment, args['command'], output_id)
        self.save(f'command-{output_id:03d}-launch.json', self.commands[output_id])
        return await self.poll_command(output_id, environment, 20)

    async def poll_command(self, output_id, environment, wait_seconds):
        if self.outputs[output_id] is not None:
            return {'state': 'complete', 'command_id': output_id, **self.page(output_id, 0)}
        info = self.commands[output_id]
        state = await durable_shell.poll(environment, info, wait_seconds)
        self.save(f'command-{output_id:03d}-status.json', state)
        if state['state'] not in ('complete', 'failed'):
            return {'state': 'running', 'command_id': output_id,
                    **{k: state[k] for k in ('stdout_path', 'stderr_path', 'stdout_bytes', 'stderr_bytes')},
                    'instruction': 'Poll this command_id; do not restart it. Complete output will be paginated after completion.'}
        record = {'command': info['command'], 'stdout': state['stdout'],
                  'stderr': state['stderr'], 'return_code': state.get('return_code'),
                  'state': state['state'], 'command_id': output_id,
                  'started_epoch': state.get('started_epoch'), 'ended_epoch': state.get('ended_epoch')}
        if state['state'] == 'failed':
            record['runner_error'] = state.get('error_type', 'unknown')
        self.save(f'command-{output_id:03d}.json', record)
        self.outputs[output_id] = json.dumps(record, ensure_ascii=False)
        return {'state': state['state'], 'command_id': output_id, **self.page(output_id, 0)}

    async def run(self, instruction, environment, context):
        model = self.model_name.split('/')[-1]
        if model not in (*CLAUDE_RATES, *OPENAI_RATES, 'gemini-3.1-pro-preview'):
            raise ValueError('This harness requires the declared paper baseline models')
        claude = model.startswith('claude')
        gemini = model.startswith('gemini')
        key = credentials('ANTHROPIC_API_KEY' if claude else 'GEMINI_API_KEY' if gemini else 'OPENAI_API_KEY')
        if gemini:
            messages = [{'role': 'user', 'parts': [{'text': instruction}]}]
        else:
            messages = [{'role': 'user', 'content': instruction}]
        retries = []
        context.n_input_tokens = context.n_output_tokens = context.n_cache_tokens = 0
        context.cost_usd = 0.0
        context.metadata = {'completed': False, 'truncation_events': 0, 'api_calls': 0,
                            'model': model, 'max_cost_usd': self.max_cost,
                            'reasoning': 'adaptive-high' if claude else 'thinking-level-high' if gemini else 'high',
                            'provider_retries': 0,
                            'prompt_cache_policy': 'automatic-5m' if claude and self.prompt_cache else 'provider-default'}
        self.save('harness-config.json', {'system': self.system, 'tools': TOOLS,
                   'instruction': instruction, 'model': model, 'max_calls': self.max_calls,
                   'max_cost_usd': self.max_cost, 'max_output_tokens': self.max_output,
                   'harness_version': self.version(), 'shell_output_policy': 'durable commands; full stdout/stderr retained'})
        for turn in range(self.max_calls):
            if claude:
                payload = {'model': model, 'max_tokens': self.max_output, 'system': self.system,
                           'messages': messages, 'thinking': {'type': 'adaptive'},
                           'output_config': {'effort': 'high'},
                           'tools': [{'name': t['name'], 'description': t['description'],
                                      'input_schema': t['parameters']} for t in TOOLS]}
                url = 'https://api.anthropic.com/v1/messages'
                headers = {'x-api-key': key, 'anthropic-version': '2023-06-01'}
                if self.prompt_cache:
                    payload['cache_control'] = {'type': 'ephemeral'}
            elif gemini:
                # Function-call parts keep their thought signatures because the
                # model's content is appended to the history verbatim.
                payload = {'systemInstruction': {'parts': [{'text': self.system}]}, 'contents': messages,
                           'tools': [{'functionDeclarations': [
                               {'name': t['name'], 'description': t['description'],
                                'parametersJsonSchema': t['parameters']} for t in TOOLS]}],
                           'generationConfig': {'maxOutputTokens': self.max_output,
                                                'thinkingConfig': {'thinkingLevel': 'high'}}}
                url = f'https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent'
                headers = {'x-goog-api-key': key}
            else:
                payload = {'model': model, 'max_output_tokens': self.max_output, 'instructions': self.system,
                           'input': messages, 'reasoning': {'effort': 'high'}, 'store': False,
                           'tools': [{'type': 'function', 'strict': True, **t} for t in TOOLS]}
                url = 'https://api.openai.com/v1/responses'
                headers = {'Authorization': 'Bearer ' + key}
            # Conservative pre-call bound: UTF-8 bytes upper-bound text tokens;
            # count the serialized request plus overhead and maximum output.
            upper_input = len(json.dumps(payload, ensure_ascii=False).encode()) + 4096
            if claude:
                rate = CLAUDE_RATES[model]
                input_rate, output_rate = (rate[2] if self.prompt_cache else rate[0]), rate[3]
            elif gemini:
                # Long-context tier applies to prompts above 200k tokens.
                input_rate, output_rate = (2, 12) if upper_input <= 200000 else (4, 18)
            else:
                rate = OPENAI_RATES[model]
                long_rate = rate[4] if rate[3] and upper_input > rate[3] else None
                input_rate, output_rate = (long_rate[0], long_rate[2]) if long_rate else (rate[0], rate[2])
            upper_cost = (upper_input * input_rate + self.max_output * output_rate) / 1e6
            if context.cost_usd + upper_cost > self.max_cost:
                self.stop_at_resource_limit(context, 'monetary_cap_before_next_request')
                return
            self.save(f'request-{turn:03d}.json', payload)
            # Non-streaming: the socket is idle until the whole response is generated.
            response = await asyncio.to_thread(post, url, headers, payload, retries,
                                               max(240, self.max_output // 30))
            context.metadata['provider_retries'] = len(retries)
            if retries:
                self.save('provider-retries.json', retries)
            self.save(f'response-{turn:03d}.json', response)
            usage = response['usageMetadata'] if gemini else response['usage']
            n_messages = len(messages)
            if claude:
                cached = usage.get('cache_read_input_tokens', 0)
                write = usage.get('cache_creation_input_tokens', 0)
                inp, out = usage['input_tokens'], usage['output_tokens']
                rate = CLAUDE_RATES[model]
                cost = (inp * rate[0] + cached * rate[1] + write * rate[2] + out * rate[3]) / 1e6
                inp += cached + write
                truncated = response['stop_reason'] == 'max_tokens'
                content = response['content']
                calls = [c for c in content if c['type'] == 'tool_use']
                normal_end = response['stop_reason'] == 'end_turn'
                messages.append({'role': 'assistant', 'content': content})
            elif gemini:
                inp = usage.get('promptTokenCount', 0)
                cached = usage.get('cachedContentTokenCount', 0)
                # Billed output includes thinking tokens.
                out = usage.get('candidatesTokenCount', 0) + usage.get('thoughtsTokenCount', 0)
                long_context = inp > 200000
                cost = ((inp-cached)*(4 if long_context else 2) + cached*(.4 if long_context else .2)
                        + out*(18 if long_context else 12)) / 1e6
                candidate = (response.get('candidates') or [{}])[0]
                truncated = candidate.get('finishReason') == 'MAX_TOKENS'
                content = candidate.get('content') or {}
                calls = [part['functionCall'] for part in content.get('parts', []) if 'functionCall' in part]
                normal_end = candidate.get('finishReason') == 'STOP'
                messages.append({'role': 'model', 'parts': content.get('parts', [])})
            else:
                inp, out = usage['input_tokens'], usage['output_tokens']
                cached = usage.get('input_tokens_details', {}).get('cached_tokens', 0)
                rate = OPENAI_RATES[model]
                if rate[3] and inp > rate[3]:
                    rate = rate[4]
                cost = ((inp-cached)*rate[0] + cached*rate[1] + out*rate[2]) / 1e6
                truncated = response.get('status') == 'incomplete'
                content = response.get('output', [])
                calls = [c for c in content if c['type'] == 'function_call']
                normal_end = response.get('status') == 'completed'
                messages.extend(content)
            context.n_input_tokens += inp
            context.n_output_tokens += out
            context.n_cache_tokens += cached
            context.cost_usd += cost
            context.metadata['api_calls'] += 1
            context.metadata['truncation_events'] += int(truncated)
            self.save('usage.json', context.model_dump())
            if truncated and context.metadata.get('truncation_retries', 0) < self.truncation_retries:
                context.metadata['truncation_retries'] = context.metadata.get('truncation_retries', 0) + 1
                del messages[n_messages:]
                self.save('usage.json', context.model_dump())
                continue
            if truncated:
                raise RuntimeError('Incomplete provider output: smoke fails, full response retained')
            if not calls:
                if not normal_end:
                    raise RuntimeError('Provider stopped without a normal completed turn')
                for command_id in self.commands:
                    if self.outputs[command_id] is None:
                        await self.poll_command(command_id, environment, 0)
                if any(output is None for output in self.outputs):
                    self.stop_at_resource_limit(context, 'provider_finished_with_active_shell_commands')
                    return
                context.metadata['completed'] = True
                context.metadata['termination_reason'] = 'normal_provider_completion'
                self.save('usage.json', context.model_dump())
                return
            results = []
            for call in calls:
                name = call['name']
                args = call['input'] if claude else call.get('args') or {} if gemini else json.loads(call['arguments'])
                result = json.dumps(await self.execute(name, args, environment), ensure_ascii=False)
                if claude:
                    results.append({'type': 'tool_result', 'tool_use_id': call['id'], 'content': result})
                elif gemini:
                    # The same serialized tool result string the other providers receive.
                    response_part = {'name': name, 'response': {'result': result}}
                    if 'id' in call:
                        response_part['id'] = call['id']
                    results.append({'functionResponse': response_part})
                else:
                    messages.append({'type': 'function_call_output', 'call_id': call['call_id'], 'output': result})
            if claude:
                messages.append({'role': 'user', 'content': results})
            elif gemini:
                messages.append({'role': 'user', 'parts': results})
        self.stop_at_resource_limit(context, 'maximum_api_calls')
