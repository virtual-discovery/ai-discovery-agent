"""Pass rates of the correct method and of the wrong paths under other budgeted purchase plans.

    cd benchmarks/construction/hzo-vacancy-local-environment && python3 -B calibrate.py

Re-uses make_task.py's payload builders and task_lab.grade's arithmetic (via make_task.evaluate and the
thresholds in hidden.json). Each design draws the purchased sites at random (class-targeted like the oracle, or
uniformly at random among the purchasable sites, with n sites per composition) and 10 random CASTEP-LDA sites.
"""
import json
import math
import random
import numpy as np
import model as M
import make_task as T

HIDDEN = json.loads((T.TASK/'environment/lab/private/hidden.json').read_text())
AVAIL = {x: [s for s in sorted(M.EF[x]) if s not in T.HELD[x]] for x in ('0.25', '0.75')}
LDA = {r['site']: r for r in M.POOL['castep']['LDA']}


def passes(parts):
    e = T.evaluate(parts, HIDDEN)
    t = HIDDEN['thresholds']
    ok = {}
    for x in ('x0.25', 'x0.75'):
        ok['pred_' + x] = e[x]['rmse'] <= t['pred_rmse_max_eV'] and abs(e[x]['bias']) <= t['pred_bias_max_eV']
        ok['dist_' + x] = (abs(e[x]['mean_nZr'] - HIDDEN['reference'][x]['mean_nZr']) <= t['mean_nZr_tol'][x]
                           and abs(e[x]['share_4Zr'] - HIDDEN['reference'][x]['share_4Zr']) <= t['share_4Zr_tol'][x])
        ok['sf_' + x] = abs(e[x]['dlog10_sf']) <= t['log10_site_fraction_tol'][x]
    ok['lda'] = abs(e['lda']['dD']) <= t['lda_offset_tol_eV'] and abs(e['lda']['db']) <= t['lda_slope_tol_eV']
    return ok


def design(seed, mode, n):
    r = random.Random(seed)
    bought = {}
    for x in ('0.25', '0.75'):
        if mode == 'targeted':
            chosen = []
            for cls, k in T.PLAN[x]:
                cand = [s for s in AVAIL[x] if M.ENV[x][s] == cls]
                chosen += r.sample(cand, k)
            bought[x] = chosen
        else:
            bought[x] = r.sample(AVAIL[x], n)
    lda = r.sample(sorted(LDA), 10)
    return bought, lda


def correct_payload(bought, lda, kind='local_x'):
    train = T.rows_for('0.50', sorted(M.EF['0.50'])) + T.rows_for('0.25', bought['0.25']) + T.rows_for('0.75', bought['0.75'])
    coef = M.fit(train, kind)
    fn = lambda x, cls: M.predict(coef, float(x), cls, kind)
    resid = np.array([e - M.predict(coef, x, cls, kind) for x, cls, e in train])
    jensen = float(np.mean(resid ** 2)) / (2 * M.KT ** 2)
    parts = T.payload_parts(fn, T.pops_from_class_energy(fn, extra=jensen))
    parts['lda_common_slope_model'] = T.lda_field(*T.common_slope(T.lda_items([LDA[s] for s in lda])))
    return parts


def blind_payload(bought, lda):
    return correct_payload(bought, lda, kind='local')


def main():
    for label, mode, n, builder in [('correct, class-targeted 9+9', 'targeted', 9, correct_payload),
                                    ('correct, random 9+9', 'random', 9, correct_payload),
                                    ('correct, random 6+6', 'random', 6, correct_payload),
                                    ('correct, random 12+12', 'random', 12, correct_payload),
                                    ('composition-blind, class-targeted', 'targeted', 9, blind_payload),
                                    ('composition-blind, random 12+12', 'random', 12, blind_payload)]:
        rows = []
        for seed in range(400):
            b, l = design(seed, mode, n)
            try:
                rows.append(passes(builder(b, l)))
            except np.linalg.LinAlgError:
                continue
        keys = rows[0].keys()
        rate = {k: round(sum(r[k] for r in rows) / len(rows), 3) for k in keys}
        allp = round(sum(all(r.values()) for r in rows) / len(rows), 3)
        print('%-36s all-pass %.3f  %s' % (label, allp, rate))


def lda_section():
    t, ref = HIDDEN['thresholds'], HIDDEN['reference']['lda']
    ok = lambda D, b: abs(D - ref['offset_O3_eV']) <= t['lda_offset_tol_eV'] and abs(b - ref['slope_per_Zr_eV']) <= t['lda_slope_tol_eV']
    for n in (6, 8, 10, 12, 14):
        good = bad = 0
        for seed in range(1000):
            sites = random.Random(seed).sample(sorted(LDA), n)
            rows = [LDA[s] for s in sites]
            try:
                good += ok(*T.common_slope(T.lda_items(rows)))
                bad += ok(*T.common_slope(T.lda_items(rows, 'vasp_cell')))
            except np.linalg.LinAlgError:
                pass
        print('LDA %2d random sites: correct labels pass %.3f | VASP-cell environments pass %.3f' % (n, good / 1000, bad / 1000))


if __name__ == '__main__':
    main()
    lda_section()
