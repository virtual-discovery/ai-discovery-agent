"""Check eliashberg.py against the authors' archived Tc tables (fig2e/tc_*.dat: McMillan, Allen-Dynes and four
Eliashberg variants for the FA spectra at mu* = 0-0.20). Writes solver_validation.json.

    python3 -B validate_solver.py      (run in the runtime image; about 10 minutes)
"""
import json
from pathlib import Path
import numpy as np
import eliashberg as E

HERE = Path(__file__).resolve().parent
POOL = json.loads((HERE/'pool.json').read_text())


def main():
    out, worst = [], 0.0
    for tag, table in sorted(POOL['tc_tables_fig2e'].items()):
        mu = int(tag.split('_')[1]) / 100
        for r in table['rows']:
            soc, x = int(r[0]), '%.2f' % r[1]
            phase = '1x1H-SOC' if soc else '1x1H'
            rows = np.array(POOL['spectra'][f'{phase}/{x}/FA']['rows'])
            s = E.Spectrum(rows[:, 0], rows[:, 1])
            d = POOL['dos'][f'{phase}/{x}']
            dos = E.DOS(d['energy_eV'], d['total'], float(x))
            ein = E.Einstein(s.lam, s.wlog)
            got = {'MM': E.allen_dynes(s.lam, s.wlog, mu, correct=False),
                   'AD': E.allen_dynes(s.lam, s.wlog, mu, s.w2),
                   'CE': E.critical_temperature(ein, mu)[0], 'CA': E.critical_temperature(s, mu)[0],
                   'DE': E.critical_temperature(ein, mu, dos)[0], 'DA': E.critical_temperature(s, mu, dos)[0]}
            ref = dict(zip(['MM', 'AD', 'CE', 'CA', 'DE', 'DA'], r[2:8]))
            dev = {k: got[k] / ref[k] - 1 for k in ref}
            worst = max(worst, max(abs(v) for v in dev.values()))
            out.append({'mu_star': mu, 'soc': bool(soc), 'x': float(x), 'archived': ref, 'this_solver': got,
                        'relative_deviation': dev})
            print(tag, phase, x, ' '.join('%s %+.4f%%' % (k, 100 * v) for k, v in dev.items()), flush=True)
    (HERE/'solver_validation.json').write_text(json.dumps({'n_values': 6 * len(out), 'max_abs_relative_deviation': worst,
                                                           'rows': out}, indent=1) + '\n')
    print('max |deviation| = %.4f %%' % (100 * worst))


if __name__ == '__main__':
    main()
