"""Host orchestration of resumable task-container shell commands; no provider calls."""
import base64,json,shlex,uuid
RUNNER = r'''
from pathlib import Path
import json,subprocess,time,os,sys,traceback
root=Path(sys.argv[1]);spec=json.loads((root/'spec.json').read_text())
def save(state):
 tmp=root/'state.tmp';tmp.write_text(json.dumps(state));tmp.replace(root/'state.json')
state={'state':'running','runner_pid':os.getpid(),'started_epoch':time.time()};save(state)
try:
 with (root/'stdout.log').open('wb') as out,(root/'stderr.log').open('wb') as err:
  p=subprocess.Popen(['bash','-lc',spec['command']],cwd='/app',stdout=out,stderr=err)
  state['command_pid']=p.pid;save(state);code=p.wait()
 state.update(state='complete',return_code=code,ended_epoch=time.time());save(state)
except BaseException as exc:
 with (root/'stderr.log').open('a') as f:traceback.print_exc(file=f)
 state.update(state='failed',error_type=type(exc).__name__,ended_epoch=time.time());save(state)
'''
async def launch(environment,command,command_id):
 root='/app/.benchmark-shell/'+uuid.uuid4().hex
 code='''from pathlib import Path
import json,subprocess,base64,sys
root=Path(ROOT);root.mkdir(parents=True,exist_ok=False)
(root/'spec.json').write_text(json.dumps({'command':COMMAND}))
(root/'runner.py').write_bytes(base64.b64decode(RUNNER64))
with (root/'runner-launch.log').open('wb') as out:
 p=subprocess.Popen([sys.executable,str(root/'runner.py'),str(root)],stdin=subprocess.DEVNULL,stdout=out,stderr=out,start_new_session=True)
print(json.dumps({'root':str(root),'runner_pid':p.pid}))
'''
 code='ROOT='+repr(root)+'\nCOMMAND='+repr(command)+'\nRUNNER64='+repr(base64.b64encode(RUNNER.encode()).decode())+'\n'+code
 result=await environment.exec('python3 -c '+shlex.quote(code),cwd='/app',timeout_sec=30)
 if result.return_code:raise RuntimeError('Could not launch durable shell: '+(result.stderr or '')[-500:])
 info=json.loads(result.stdout);return {**info,'command':command,'command_id':command_id}

async def poll(environment,info,wait_seconds):
 code='''from pathlib import Path
import json,time,os
root=Path(ROOT);deadline=time.monotonic()+WAIT
while True:
 try:state=json.loads((root/'state.json').read_text())
 except FileNotFoundError:state={'state':'starting'}
 if state['state'] in ('complete','failed'):break
 if time.monotonic()>=deadline:break
 time.sleep(.2)
state['stdout_path']=str(root/'stdout.log');state['stderr_path']=str(root/'stderr.log')
for stream in ('stdout','stderr'):
 p=root/(stream+'.log');state[stream+'_bytes']=p.stat().st_size if p.exists() else 0
if state['state'] in ('complete','failed'):
 state['stdout']=(root/'stdout.log').read_text(errors='replace') if (root/'stdout.log').exists() else ''
 state['stderr']=(root/'stderr.log').read_text(errors='replace') if (root/'stderr.log').exists() else ''
print(json.dumps(state))
'''
 code='ROOT='+repr(info['root'])+'\nWAIT='+repr(wait_seconds)+'\n'+code
 result=await environment.exec('python3 -c '+shlex.quote(code),cwd='/app',timeout_sec=wait_seconds+20)
 if result.return_code:raise RuntimeError('Durable shell status query failed; command state unknown')
 return json.loads(result.stdout)
