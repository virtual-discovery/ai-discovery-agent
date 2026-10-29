"""Reference SSMF evaluations through the lab path (kag_ssmf.lab_ssmf). Appends to /data/refs.json.

    bash run_in_image.sh make_refs.py MODEL:N:U:JU[:mesh] ...
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, '.')
import kag_ssmf as K
FOLDER = {'TiSb': 'CsTi3Sb5', 'VSb': 'CsV3Sb5', 'CrSb': 'CsCr3Sb5', 'TiBi': 'CsTi3Bi5', 'CrBi': 'CsCr3Bi5',
          'CrBi_9.2': 'CsCr3Bi5_9.2'}
out = Path(sys.argv[1])
for job in sys.argv[2:]:
    parts = job.split(':')
    key, N, U, JU = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
    m = K.load_model(f'/data/models/{FOLDER[key]}.npz')
    s = K.SSMF(m, mesh=tuple(int(x) for x in parts[4].split(',')) if len(parts) > 4 else K.LAB_MESH, kT=K.LAB_KT)
    t = time.time()
    r = K.lab_ssmf(m, N, U, JU, solver=s)
    r['seconds'] = round(time.time() - t, 1)
    data = json.loads(out.read_text()) if out.exists() else {}
    data[job] = r
    out.write_text(json.dumps(data, indent=1))
    print(job, r['metallic'], np.round(r.get('Z', []), 4), r['seconds'], flush=True)
