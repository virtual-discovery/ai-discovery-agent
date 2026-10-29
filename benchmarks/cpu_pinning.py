"""Give each Quantum ESPRESSO lab container its own physical cores.

The lab services launch `mpirun -np N pw.x`. OpenMPI binds rank i to core i
inside every container, so concurrent fresh/live episodes all compete for host
cores 0..N-1 while the rest of the machine idles. Slower SCFs then cost
agents extra polling calls against their call and dollar budgets. This host
daemon enforces each task's declared `cpus` with dedicated physical cores
(including hyperthread siblings) through `docker update --cpuset-cpus`. A
cpuset change also re-spreads ranks that are already bound. Task packages are
unchanged. Run one instance per host during any campaign with QE tasks:

    python -m benchmarks.cpu_pinning
"""
import fcntl
import json
from pathlib import Path
import re
import subprocess
import sys
import time
from benchmarks.verify import ROOT

TASKS = ROOT/'benchmarks/tasks'
HEARTBEAT = ROOT/'runs/.cpu-pinning.heartbeat'
LOG = ROOT/'runs/cpu-pinning.log'
RESERVED_CORES = 2   # physical cores 0-1 of package 0 stay with the host and Harbor


def qe_tasks():
    tasks = {}
    for folder in [f for f in TASKS.iterdir() if f.is_dir()]:
        sources = ' '.join(p.read_text() for p in (folder/'environment/lab').glob('*.py')) if (folder/'environment/lab').is_dir() else ''
        if 'pw.x' in sources:
            cpus = re.search(r'^cpus\s*=\s*(\d+)', (folder/'task.toml').read_text(), re.M)
            tasks[folder.name] = int(cpus.group(1)) if cpus else 2
    return tasks


def physical_cores():
    cores = {}
    for cpu in Path('/sys/devices/system/cpu').glob('cpu[0-9]*'):
        topology = cpu/'topology'
        key = (int((topology/'physical_package_id').read_text()), int((topology/'core_id').read_text()))
        cores.setdefault(key, set()).add(int(cpu.name[3:]))
    ordered = sorted(cores)
    return {key: sorted(cores[key]) for key in ordered[RESERVED_CORES:]}


def log(entry):
    with LOG.open('a') as stream:
        stream.write(json.dumps({'epoch': time.time(), **entry}) + '\n')


def main():
    LOG.parent.mkdir(exist_ok=True)
    lock = (ROOT/'runs/.cpu-pinning.lock').open('w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit('Another cpu_pinning daemon is running')
    tasks, cores = qe_tasks(), physical_cores()
    assigned = {}   # container id -> list of physical-core keys
    log({'event': 'start', 'qe_tasks': tasks, 'cores': len(cores)})
    scanned = time.time()
    while True:
        if time.time() - scanned > 30:   # tasks built while the daemon runs
            tasks, scanned = qe_tasks(), time.time()
        listing = subprocess.run(['docker', 'ps', '--format', '{{.ID}} {{.Names}}'],
                                 capture_output=True, text=True).stdout.split('\n')
        running = dict(line.split(' ', 1) for line in listing if ' ' in line)
        for gone in set(assigned) - set(running):
            log({'event': 'release', 'container': gone, 'cores': assigned.pop(gone)})
        for cid, name in running.items():
            task = name.split('__')[0]
            if cid in assigned or not name.endswith('-lab-1') or task not in tasks:
                continue
            used = {key for keys in assigned.values() for key in keys}
            free = [key for key in cores if key not in used]
            # Prefer a single socket for one container's ranks.
            by_socket = [[k for k in free if k[0] == socket] for socket in sorted({k[0] for k in free})]
            pick = next((group[:tasks[task]] for group in by_socket if len(group) >= tasks[task]), free[:tasks[task]])
            if len(pick) < tasks[task]:
                log({'event': 'no_free_cores', 'container': name})
                continue
            cpuset = ','.join(str(c) for key in pick for c in cores[key])
            result = subprocess.run(['docker', 'update', '--cpuset-cpus', cpuset, cid], capture_output=True, text=True)
            if result.returncode == 0:
                assigned[cid] = pick
                log({'event': 'pin', 'container': name, 'cpuset': cpuset})
            else:
                log({'event': 'pin_failed', 'container': name, 'error': result.stderr[-300:]})
        HEARTBEAT.write_text(json.dumps({'epoch': time.time(), 'assigned': {running.get(c, c): k for c, k in assigned.items()}}))
        time.sleep(1)


if __name__ == '__main__':
    main()
