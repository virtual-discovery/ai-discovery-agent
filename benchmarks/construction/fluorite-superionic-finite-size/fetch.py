"""Refetch the pinned source files for fluorite-superionic-finite-size by HTTP range reads.

    python3 fetch.py            # -> ~/.cache/ai-theorist/fluorite-superionic-finite-size/raw/

Source: F. Grasselli, Materials Cloud Archive record yzhmw-vje16 (doi:10.24435/materialscloud:jy-tw), CC BY 4.0.
Only the diffusion (NVT, MSD) folders are read, never the *_heat folders (those feed the sibling task
superionic-heat-transport-size-onsager). Members are read out of the remote zips without downloading the
archives; each member is checked against the zip CRC32, and its SHA256 is compared with members.json when present.
nvt_log files (6.6 MB text each) are reduced on the fly to float64 .npy arrays of the columns
(Time [ps], TotEng, PotEng [eV], Temp [K], Press [bar]); the SHA256 recorded is that of the original text member.
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
RAW = Path(os.path.expanduser('~/.cache/ai-theorist/fluorite-superionic-finite-size/raw'))
ARCHIVES = {'CaF2.zip': '7a88bad620f87627aa4acb18266254d1', 'PbF2.zip': '7c99d63a2e0e0662a2aa634496405e39',
            'UO2.zip': 'ba43b29a5c7308bcca86b29041c559a6', 'alpha-AgI.zip': '6532ccd122eaf453975676e0e8acbbc0'}
SETS = {'CaF2.zip': ['CaF2/NVT_L_T/'], 'PbF2.zip': ['PbF2/NVT_L_T_6p056/', 'PbF2/NVT_L_T_6p056_DW/'],
        'UO2.zip': ['UO2/NVT_L_T/'], 'alpha-AgI.zip': ['alpha-AgI/NVT_L_T_5.37ang/']}
KEEP = ('msd_CM_B4.dat', 'msd_self_B4.dat', 'nvt_log', '.in', 'run_all_T_and_L.sh', 'run_msd_analisi.sh',
        'run_extract_log.sh')


def wanted(name, prefixes):
    return any(name.startswith(p) for p in prefixes) and name.endswith(KEEP) and not name.endswith('/')


def get(zipname, entry, pinned):
    out = RAW/entry['name']
    target = out.with_suffix('.npy') if out.name == 'nvt_log' else out
    if target.exists() and entry['name'] in pinned:
        return entry['name'], pinned[entry['name']]
    data = RZ.extract(zipname, entry)
    digest = hashlib.sha256(data).hexdigest()
    if entry['name'] in pinned and pinned[entry['name']] != digest:
        raise RuntimeError(f"SHA256 mismatch for {entry['name']}")
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.name == 'nvt_log':
        lines = data.decode().split('\n')
        assert lines[0].split()[:7] == ['Step', 'Time', 'TotEng', 'PotEng', 'Temp', 'Press', 'Lx'], lines[0]
        arr = np.array([[float(x) for x in l.split()[1:6]] for l in lines[1:] if l.strip()])
        np.save(target, arr)
    else:
        out.write_bytes(data)
    return entry['name'], digest


def main():
    pinned_path = HERE/'members.json'
    pinned = json.loads(pinned_path.read_text())['sha256'] if pinned_path.exists() else {}
    record = {}
    for zipname, prefixes in SETS.items():
        entries = [e for e in RZ.central_directory(zipname) if wanted(e['name'], prefixes)]
        print(zipname, len(entries), 'members', flush=True)
        with ThreadPoolExecutor(4) as pool:
            for name, digest in pool.map(lambda e: get(zipname, e, pinned), entries):
                record[name] = digest
    if not pinned:
        pinned_path.write_text(json.dumps({'archive_md5': ARCHIVES, 'record': RZ.RECORD,
                                           'sha256': dict(sorted(record.items()))}, indent=0) + '\n')
    missing = set(pinned) - set(record)
    if missing:
        sys.exit(f'members missing from the archive: {sorted(missing)[:5]}')
    print(len(record), 'members verified ->', RAW)


if __name__ == '__main__':
    main()
