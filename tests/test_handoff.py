import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from benchmarks.run import select, build_command, completed_job, cost_in_runs
from benchmarks.tools.release import task_rows
from benchmarks.api_agent import credentials, PaperBaselineAgent
from harbor.models.agent.context import AgentContext


class HandoffTests(unittest.TestCase):
    def test_v2_selection_and_limits(self):
        rows = task_rows()
        self.assertEqual(len(rows), 60)
        self.assertEqual(len(select(rows, 'all')), 60)
        with self.assertRaises(ValueError):
            select(rows, 'all', tasks=['not-a-task'])
        _, command = build_command({'task_id': 'heis-configdesign-fe', 'start': 'released'}, 'astra', Path('/tmp/r'), 1, 20)
        for flag in ('max_calls=200', 'max_cost_usd=20', 'agent_qe=false', 'max_output_tokens=32000'):
            self.assertIn(flag, command)

    def test_failed_terminal_is_not_retried_and_incomplete_cost_is_counted(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            job = root/'job'
            usage = job/'trial/agent/usage.json'
            usage.parent.mkdir(parents=True)
            usage.write_text(json.dumps({'cost_usd': .7, 'metadata': {'completed': False}}))
            self.assertAlmostEqual(cost_in_runs(root), .7)
            self.assertFalse(completed_job(job))
            (job/'result.json').write_text(json.dumps({'finished_at': '2026-10-09T00:00:00Z', 'reward': 0}))
            self.assertTrue(completed_job(job))

    def test_host_credentials_need_no_workstation_file(self):
        with patch.dict('os.environ', {'OPENAI_API_KEY': 'unit-test-placeholder'}):
            self.assertEqual(credentials('OPENAI_API_KEY'), 'unit-test-placeholder')
        with tempfile.TemporaryDirectory() as temp, patch('benchmarks.api_agent.ROOT', Path(temp)), patch.dict('os.environ', {}, clear=True):
            with self.assertRaisesRegex(RuntimeError, 'Set OPENAI_API_KEY'):
                credentials('OPENAI_API_KEY')

    def test_gemini_history_keeps_signatures_and_costs_thinking(self):
        responses = [
            {'candidates': [{'finishReason': 'STOP', 'content': {'role': 'model', 'parts': [
                {'functionCall': {'name': 'shell', 'args': {'command': 'ls'}, 'id': 'call_1'}, 'thoughtSignature': 'sig-a'}]}}],
             'usageMetadata': {'promptTokenCount': 1000, 'cachedContentTokenCount': 400, 'candidatesTokenCount': 10, 'thoughtsTokenCount': 90}},
            {'candidates': [{'finishReason': 'STOP', 'content': {'role': 'model', 'parts': [{'text': 'done', 'thoughtSignature': 'sig-b'}]}}],
             'usageMetadata': {'promptTokenCount': 1200, 'candidatesTokenCount': 5}}]
        payloads = []

        def fake_post(url, headers, payload, retries=None, timeout=None):
            payloads.append(json.loads(json.dumps(payload)))
            self.assertIn('gemini-3.1-pro-preview:generateContent', url)
            return responses[len(payloads)-1]

        async def fake_execute(name, args, environment):
            return {'state': 'complete', 'output': 'listing'}
        with tempfile.TemporaryDirectory() as temp, patch('benchmarks.api_agent.post', fake_post), \
                patch.dict('os.environ', {'GEMINI_API_KEY': 'unit-test-placeholder'}):
            agent = PaperBaselineAgent(logs_dir=Path(temp), model_name='google/gemini-3.1-pro-preview')
            agent.execute = fake_execute
            context = AgentContext()
            import asyncio
            asyncio.run(agent.run('task', None, context))
        self.assertTrue(context.metadata['completed'])
        second = payloads[1]['contents']
        self.assertEqual(second[1]['parts'][0]['thoughtSignature'], 'sig-a')
        reply = second[2]['parts'][0]['functionResponse']
        self.assertEqual((reply['id'], reply['name']), ('call_1', 'shell'))
        self.assertEqual(json.loads(reply['response']['result'])['output'], 'listing')
        self.assertEqual(payloads[0]['generationConfig']['thinkingConfig'], {'thinkingLevel': 'high'})
        expected = (600*2 + 400*.2 + 100*12 + 1200*2 + 5*12) / 1e6
        self.assertAlmostEqual(context.cost_usd, expected)
        self.assertEqual(context.n_output_tokens, 105)


if __name__ == '__main__':
    unittest.main()
