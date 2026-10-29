"""Refetch the pinned heat-transport members for superionic-heat-transport-size-onsager by HTTP range reads.

    taskset -c 6,7,22,23 python3 -B fetch.py     # -> ~/.cache/ai-theorist/superionic-heat-transport-size-onsager/raw/

Source: F. Grasselli, Materials Cloud Archive record yzhmw-vje16 (doi:10.24435/materialscloud:jy-tw), CC BY 4.0.
Only the heat-transport folders are read (never the NVT diffusion folders used by the sibling task
fluorite-superionic-finite-size):
  UO2_heat.zip: UO2_heat/NVE_L_T/, UO2_heat/NVT_L_T_heat/
  CaF2.zip:     CaF2/NVE_L_T_heat/, CaF2/NVT_L_T_heat/
  PbF2.zip:     PbF2/NVE_L_T_6p056_heat/, PbF2/NVT_L_T_6p056_heat/
Members are read out of the remote zips without downloading the archives; each member is checked against
the zip CRC32 and its SHA256 is compared with members.json when present. Each thermo log (nve_log/nvt_log,
24 MB text, 100001 rows every 8 fs) is reduced on the fly to a float32 .npy array of the columns
  Time [ps], TotEng, PotEng [eV], Temp [K], Press [bar], Lx [A], J_x, J_y, J_z [eV A/ps], V_x, V_y, V_z [A/ps]
where J is LAMMPS compute heat/flux (extensive, J*V) and V is the sum of the anion velocities (O or F).
The six stress columns are dropped. The SHA256 recorded is that of the original text member.
"""
import hashlib
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import remotezip as RZ

HERE = Path(__file__).resolve().parent
RAW = Path(os.path.expanduser('~/.cache/ai-theorist/superionic-heat-transport-size-onsager/raw'))
ARCHIVES = {'UO2_heat.zip': 'de844afbecd3fd5e72812ccb341644ad', 'CaF2.zip': '7a88bad620f87627aa4acb18266254d1',
            'PbF2.zip': '7c99d63a2e0e0662a2aa634496405e39'}
SETS = {'UO2_heat.zip': ['UO2_heat/NVE_L_T/', 'UO2_heat/NVT_L_T_heat/'],
        'CaF2.zip': ['CaF2/NVE_L_T_heat/', 'CaF2/NVT_L_T_heat/'],
        'PbF2.zip': ['PbF2/NVE_L_T_6p056_heat/', 'PbF2/NVT_L_T_6p056_heat/']}
KEEP = ('_log', '.in', '.sh')
HEADER = ['Step', 'Time', 'TotEng', 'PotEng', 'Temp', 'Press', 'Lx', 'c_myHeatFlux[1]', 'c_myHeatFlux[2]',
          'c_myHeatFlux[3]']


def wanted(name, prefixes):
    return any(name.startswith(p) for p in prefixes) and name.endswith(KEEP) and not name.endswith('/')


DONE = {}


def get(zipname, entry, pinned):
    out = RAW/entry['name']
    target = out.with_suffix('.npy') if out.name.endswith('_log') else out
    if target.exists() and entry['name'] in pinned:
        return entry['name'], pinned[entry['name']]
    if target.exists() and entry['name'] in DONE:
        return entry['name'], DONE[entry['name']]
    for attempt in range(5):
        try:
            data = RZ.extract(zipname, entry)
            break
        except Exception:
            if attempt == 4:
                raise
    digest = hashlib.sha256(data).hexdigest()
    if entry['name'] in pinned and pinned[entry['name']] != digest:
        raise RuntimeError(f"SHA256 mismatch for {entry['name']}")
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.name.endswith('_log'):
        text = data.decode()
        first, rest = text.split('\n', 1)
        cols = first.split()
        assert cols[:10] == HEADER and len(cols) == 19, cols
        rows = [l for l in rest.split('\n') if l.strip() and not l.startswith(('WARNING', 'Loop'))]
        arr = np.array(' '.join(rows).split(), dtype=np.float64).reshape(len(rows), 19)
        np.save(target, arr[:, 1:13].astype(np.float32))
    else:
        out.write_bytes(data)
    return entry['name'], digest


def main():
    pinned_path = HERE/'members.json'
    pinned = json.loads(pinned_path.read_text())['sha256'] if pinned_path.exists() else {}
    record = {}
    for f in RAW.glob('*.sha256.json'):
        DONE.update(json.loads(f.read_text()))
    only = sys.argv[1:] or list(SETS)
    for zipname in only:
        entries = [e for e in RZ.central_directory(zipname) if wanted(e['name'], SETS[zipname])]
        print(zipname, len(entries), 'members', flush=True)
        with ThreadPoolExecutor(4) as pool:
            for name, digest in pool.map(lambda e: get(zipname, e, pinned), entries):
                record[name] = digest
                DONE[name] = digest
        (RAW/f'{zipname}.sha256.json').write_text(json.dumps(record, indent=0))
    if not pinned and set(only) == set(SETS):
        pinned_path.write_text(json.dumps({'archive_md5': ARCHIVES, 'record': RZ.RECORD,
                                           'sha256': dict(sorted(record.items()))}, indent=0) + '\n')
    missing = {k for k in pinned if k.split('/')[0] + '.zip' in only or k.startswith(tuple(p for z in only for p in SETS[z]))} - set(record)
    if missing:
        sys.exit(f'members missing from the archive: {sorted(missing)[:5]}')
    print(len(record), 'members verified ->', RAW)


if __name__ == '__main__':
    main()
