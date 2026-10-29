"""Spread of the correct method over purchase plans (one SQS per composition, subsets of compositions).

    python3 -B calibrate.py
"""
import itertools, json
import numpy as np
import model as M

X7 = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
PLANS = {'7comp': X7, '5comp': [0.125, 0.25, 0.5, 0.75, 0.875], '3comp': [0.25, 0.5, 0.75], '4comp': [0.25, 0.375, 0.625, 0.75]}


def run(rng, n=400):
    rows = {fn: M.sqs_rows(fn) for fn in ('PBE', 'r2SCAN')}
    bad = {fn: M.vegard_outliers(rows[fn], fn) for fn in rows}
    valid = {fn: {} for fn in rows}
    for fn in rows:
        for r, b in zip(rows[fn], bad[fn]):
            if not b:
                valid[fn].setdefault(r['x'], []).append(r['index'])
    ac = {fn: M.strain_fit(fn, '001') for fn in rows}
    out = {}
    for plan, xs in PLANS.items():
        res = []
        for _ in range(n):
            vals = {}
            for fn in rows:
                ind = {(x, int(rng.choice(valid[fn][x]))) for x in xs}
                x, y = M.sqs_points(fn, 'clean', ind)
                a, b = M.subregular(x, y)
                vals[fn] = (a, b)
            a, b = vals['PBE']; c = (a - ac['PBE'][0], b - ac['PBE'][1])
            tc, xc = M.spinodal_Tc(*c)
            bi = M.binodal_at(*c, 700.0) or (np.nan, np.nan); sp = M.spinodal_at(*c, 700.0) or (np.nan, np.nan)
            a2, b2 = vals['r2SCAN']; tc2 = M.spinodal_Tc(a2 - ac['r2SCAN'][0], b2 - ac['r2SCAN'][1])[0]
            res.append([a / 8 * 1000, tc, xc, bi[0], bi[1], sp[0], sp[1], tc2])
        res = np.array(res)
        out[plan] = {k: [float(np.nanpercentile(res[:, i], 5)), float(np.nanpercentile(res[:, i], 50)), float(np.nanpercentile(res[:, i], 95))]
                     for i, k in enumerate(['dh05', 'TcPBE', 'xcPBE', 'bi_lo', 'bi_hi', 'sp_lo', 'sp_hi', 'Tcr2SCAN'])}
    return out


if __name__ == '__main__':
    out = run(np.random.default_rng(0))
    json.dump(out, open('calibration.json', 'w'), indent=1)
    for plan, d in out.items():
        print(plan, {k: [round(v, 3) for v in vv] for k, vv in d.items()})
