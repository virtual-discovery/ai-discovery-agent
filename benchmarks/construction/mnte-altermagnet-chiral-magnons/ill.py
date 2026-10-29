"""Parser for ILL TAS ASCII files (IN20 cryopad) used by the MnTe chiral-magnon task."""
import re
import numpy as np


def parse(path):
    lines = open(path, errors='replace').read().splitlines()
    head, i = {}, 0
    params, varia, polan = {}, {}, []
    while i < len(lines) and not lines[i].startswith('DATA_:'):
        ln = lines[i]
        if ln[:6] in ('PARAM:', 'VARIA:', 'ZEROS:'):
            for k, v in re.findall(r'(\w+)\s*=\s*([-+0-9.eE]+)', ln[6:]):
                (params if ln.startswith('PARAM') else varia if ln.startswith('VARIA') else head.setdefault('zeros', {}))[k] = float(v)
        elif ln.startswith('POLAN:'):
            polan.append(ln[6:].strip())
        elif len(ln) > 6 and ln[5] == ':':
            head[ln[:5]] = ln[6:].strip()
        i += 1
    cols = lines[i + 1].split()
    rows = [l.split() for l in lines[i + 2:] if l.strip() and not l.startswith('Finished')]
    rows = [r for r in rows if len(r) == len(cols)]
    data = np.array(rows, float)
    return {'head': head, 'param': params, 'varia': varia, 'polan': polan, 'cols': cols, 'data': data}
