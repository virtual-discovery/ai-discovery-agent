"""Parse python-sscha minimization logs: per temperature section, the final RESULTS block."""
import re
import numpy as np

NUM = r'[-+]?\d+\.\d+(?:[eE][-+]?\d+)?'


def parse_log(path):
    sections, cur = [], None
    lines = [l for l in open(path, errors='ignore') if not l.startswith('Computing configuration')
             and not l.startswith('ENSEMBLE ALL')]
    i = 0
    while i < len(lines):
        l = lines[i]
        m = re.match(r'^\s+T = (\d+)\s*$', l)
        if m and i > 0 and lines[i - 1].startswith('====='):
            cur = {'T': int(m.group(1)), 'results': [], 'n_conf': None, 'steps': 0}
            sections.append(cur)
        elif cur is not None:
            m2 = re.match(r'^Force computed shape: (\d+)', l)
            if m2:
                cur['n_conf'] = int(m2.group(1))
            if l.startswith('Step ka'):
                cur['steps'] += 1
            if 'satisfy the convergence criteria' in l:
                cur['last_stop'] = 'converged'
            if 'out of the statistical sampling' in l:
                cur['last_stop'] = 'kong_liu_limit'
            if '*   RESULTS   *' in l:
                res = {'stop': cur.get('last_stop')}
                cur['last_stop'] = None
                j = i + 1
                while j < len(lines) and '*   RESULTS   *' not in lines[j] and not re.match(r'^\s+T = \d+\s*$', lines[j]):
                    s = lines[j]
                    if s.startswith('Free energy ='):
                        v = re.findall(NUM, s); res['F_meV'] = float(v[0]); res['F_err_meV'] = float(v[1])
                    elif s.startswith('Kong-Liu'):
                        res['kong_liu'] = float(re.findall(NUM, s)[0])
                    elif s.startswith('Minimization ended after'):
                        res['steps_reported'] = int(s.split()[3])
                    elif 'STRESS TENSOR' in s:
                        rows = [re.findall(NUM, lines[j + k]) for k in (1, 2, 3)]
                        res['stress_GPa'] = [[float(x) for x in r[:3]] for r in rows]
                        res['stress_err_GPa'] = [[float(x) for x in r[3:6]] for r in rows]
                    elif 'Ab initio average stress' in s:
                        rows = [re.findall(NUM, lines[j + k]) for k in (1, 2, 3)]
                        res['avg_stress_GPa'] = [[float(x) for x in r[:3]] for r in rows]
                        break
                    elif s.startswith('Struct gradient'):
                        v = re.findall(NUM, s); res['struct_grad'] = [float(v[0]), float(v[1])]
                    elif s.startswith('FC gradient'):
                        v = re.findall(NUM, s); res['fc_grad'] = [float(v[0]), float(v[1])]
                    j += 1
                cur['results'].append(res)
                i = j - 1
        i += 1
    return sections
