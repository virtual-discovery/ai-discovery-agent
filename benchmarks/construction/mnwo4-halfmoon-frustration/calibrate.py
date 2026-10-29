"""Calibration fits for mnwo4-halfmoon-frustration (run in the runtime image, ~15 min on 4 CPUs).

    python3 -B calibrate.py CACHE      (CACHE holds extracted.npz and model.json from extract.py)

Fits the purchasable archived data in several defensible ways and along the realistic wrong paths, and writes
calibration_fits.json with every graded quantity (J, D, k_LT, E_spiral - E_AF1, block predictions of the withheld
14 K map and reduced-model maps). make_task.py turns these into thresholds, ablation and robustness values.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from model import Model, truth, TABLE_I, TABLE_I_D, map14_blocks, map14_centres, ko_blocks, ko_centres, norm3, \
    half_to_full, REDUCED, CELL
from fitlib import fit_spectra, null_scan, predictions

HERE = Path(__file__).resolve().parent
CACHE = Path(sys.argv[1])
X = np.load(CACHE/'extracted.npz')
shells, J0, D0, dist = truth(CACHE/'model.json')
M = Model(shells)
ORACLE_IDX = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 88, 94, 100]
WIDE_IDX = list(range(0, 101, 3))


def spectra(idx):
    return [(X['spectra_q'][i], X['spectra_E'], X['spectra_S'][i]) for i in idx]


MAPS = {T: X['map_T%d' % T][:, :, 3] for T in (20, 30, 60)}


def pack(pred):
    return {'J': pred['J'], 'D': pred['D'], 'k': pred['k'], 'dE': pred['dE'], 'map14': pred['map14'].tolist(),
            'reduced': {m: v.tolist() for m, v in pred['reduced'].items()}}


def formfactor(qfrac):
    """Mn2+ <j0> (Brown, International Tables C) at |q| (1/A, including 2 pi)"""
    qc = 2 * np.pi * np.linalg.norm(qfrac @ np.linalg.inv(CELL).T, axis=1)
    s = qc / (4 * np.pi)
    f = 0.4220 * np.exp(-17.684 * s ** 2) + 0.5948 * np.exp(-6.005 * s ** 2) + 0.0043 * np.exp(0.609 * s ** 2) - 0.0219
    return f


def main():
    out = {}
    t0 = time.time()
    # reference (archived Hamiltonian)
    out['reference_model'] = pack(predictions(M, np.r_[J0, D0]))
    # oracle path: 20 spectra from the published (rounded) Table I start, null direction from 20 + 30 K maps
    p, c = fit_spectra(M, spectra(ORACLE_IDX), np.r_[TABLE_I, TABLE_I_D])
    out['spectra_only_fit'] = {'p': p.tolist(), 'cost': c}
    po = null_scan(M, p, {20: MAPS[20], 30: MAPS[30]})
    out['oracle'] = pack(predictions(M, po)); print('oracle', np.round(po, 4), '%.0fs' % (time.time() - t0), flush=True)
    # robustness: 34 spectra, null direction from the 20 K map only
    p2, _ = fit_spectra(M, spectra(WIDE_IDX), np.r_[TABLE_I, TABLE_I_D])
    pw = null_scan(M, p2, {20: MAPS[20]})
    out['robust_wide_spectra_20K_only'] = pack(predictions(M, pw)); print('wide', np.round(pw, 4), flush=True)
    # robustness: null direction from the 20, 30 and 60 K maps
    p3 = null_scan(M, p, {20: MAPS[20], 30: MAPS[30], 60: MAPS[60]})
    out['robust_three_maps'] = pack(predictions(M, p3)); print('three maps', np.round(p3, 4), flush=True)
    # robustness: null direction from the 11 K full-zone peak (k_LT closest to the 36x24x36 box maximum)
    kbox = np.array([-8 / 36, 0.5, 16 / 36])
    s = p[2] + p[10]; best = None
    for J3 in np.linspace(s, 0, 81):
        q = p.copy(); q[2] = J3; q[10] = s - J3
        k, _ = M.lt(q[:12], n=24)
        d = min(np.abs((k - sg * kbox) - np.round(k - sg * kbox)).max() for sg in (1, -1))
        if best is None or d < best[0]:
            best = (d, q)
    out['robust_zone_peak_null'] = pack(predictions(M, best[1])); print('zone peak', np.round(best[1], 4), flush=True)
    # ablation dp1: dispersion only, minimum-norm split of the null direction
    pm = p.copy(); pm[2] = pm[10] = (p[2] + p[10]) / 2
    out['abl_dispersion_only_min_norm'] = pack(predictions(M, pm))
    # ablation dp2: spin length 1 in the precession (energies ~ J instead of J/S)
    ps, _ = fit_spectra(M, spectra(ORACLE_IDX), np.r_[TABLE_I, TABLE_I_D] * 0.4, S=1.0)
    ps = null_scan(M, ps, {20: MAPS[20], 30: MAPS[30]})
    out['abl_spin_length_one'] = pack(predictions(M, ps)); print('S=1', np.round(ps, 4), flush=True)
    # ablation dp3: mean-field RPA (spin-length constraint dropped: lambda fixed at its Curie value 3/2)
    pr = dict(predictions(M, po, T14=14.0, lam_fixed=1.5))
    out['abl_rpa_nominal_T'] = pack(pr)
    pr17 = predictions(M, po, T14=17.0, lam_fixed=1.5)
    out['abl_rpa_lowest_valid_T'] = pack(pr17)
    # ablation dp4: isotropic model (D = 0) fitted to the same spectra
    pi_, _ = fit_spectra(M, spectra(ORACLE_IDX), np.r_[TABLE_I, 0.0], fixD=0.0)
    pi_ = null_scan(M, pi_, {20: MAPS[20], 30: MAPS[30]})
    out['abl_isotropic'] = pack(predictions(M, pi_)); print('iso', np.round(pi_, 4), flush=True)
    # ablation dp6: magnetic form factor applied to the predicted maps
    o = out['oracle']
    out['abl_form_factor'] = {'map14': (np.array(o['map14']) * formfactor(map14_centres()).reshape(27, 29) ** 2).tolist(),
                              'reduced': {m: (np.array(v) * formfactor(ko_centres()).reshape(32, 32) ** 2).tolist()
                                          for m, v in o['reduced'].items()}}
    # ablation dp7: temperature-insensitive half-moons (archived 20 K map as the 14 K map; full model for the
    # reduced models)
    full20 = norm3(M.scga(np.array(o['J']), o['D'], 20.0, ko_centres(), nbz=32)).reshape(32, 32)
    out['abl_temperature_insensitive'] = {'map14': map14_blocks(MAPS[20]).tolist(),
                                          'reduced': {m: full20.tolist() for m in REDUCED}}
    # ablation dp8: quasi-1D zigzag-chain model (J1 and J4 only, refitted to the spectra)
    pc_ = np.r_[np.array(o['J']) * 0, o['D']]; pc_[0] = o['J'][0]; pc_[3] = o['J'][3]
    chain = predictions(M, pc_)
    out['abl_chain_only'] = pack(chain)
    # non-recall route (nonrecall.py, joint DE without literature values): all seeds, and the best-objective one
    nr = {}
    for f in sorted(HERE.glob('nonrecall_seed*.json')):
        d = json.loads(f.read_text())
        nr[f.stem] = d
        out['nonrecall_' + f.stem.split('_')[1]] = pack(predictions(M, np.array(d['p'])))
    best = min(nr, key=lambda k: nr[k]['f'])
    out['robust_nonrecall_best_of_seeds'] = out['nonrecall_' + best.split('_')[1]]
    (HERE/'calibration_fits.json').write_text(json.dumps(out))
    print('done %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
