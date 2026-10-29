"""Reference solution for wurtzite-ferroelectric-switching; runs as the Harbor oracle and in check_task.

Route (no literature numbers): converged NEB runs define the barriers; the collective path is identified from
its structure; Berry-phase branches are fixed with Born-charge estimates of the polarization change between
purchased images (segments short enough that the estimate error stays far below half a quantum).
"""
import json
import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
RY = 13.605693122994
E_CHARGE = 1.602176634e-19
N_FU = 16
EVIDENCE = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def buy(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(1)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(r)
    EVIDENCE.append(r['id'])
    return r['result']


def cation_offsets(cell, species, pos):
    """z of each cation minus the mean z of its three basal N neighbours (A); the sign gives its polarity."""
    cell, pos = np.array(cell), np.array(pos)
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) @ cell
    an = [i for i, s in enumerate(species) if s == 'N']
    out = []
    for c in [i for i, s in enumerate(species) if s != 'N']:
        d = (pos[an][None] + shifts[:, None] - pos[c]).reshape(-1, 3)
        near = d[np.argsort(np.linalg.norm(d, axis=1))[:5]]
        basal = near[np.argsort(-np.linalg.norm(near[:, :2], axis=1))[:3]]
        out.append(-basal[:, 2].mean())
    return np.array(out), [s for s in species if s != 'N']


def berry_p(result):
    m = re.search(r'P =\s+(-?\d+\.\d+)\s+\(mod\s+(\d+\.\d+)\)\s+C/m\^2', result['polarization_output'])
    return float(m.group(1)), float(m.group(2))


def main():
    runs = {('AlN', 'refined_NEB2'), ('AlBN', 'NEB'), ('AlBN', 'refined_NEB'), ('AlBN', 'refined_NEB2')}
    neb = {k: buy({'kind': 'neb_run', 'system': k[0], 'run': k[1]}) for k in sorted(runs)}
    converged = {}
    for k, r in neb.items():
        thr = float(r['settings']['path_thr'].split()[0])
        converged[k] = r['convergence_message'] is not None and max(x['error_eV_per_A'] for x in r['final_table']) <= thr
    E = {k: np.array([x['energy_eV'] for x in r['final_table']]) for k, r in neb.items()}
    # Minimum-energy path = lowest-barrier converged run of each system
    mep = {}
    for s in ('AlN', 'AlBN'):
        cands = [k for k in neb if k[0] == s and converged[k]]
        mep[s] = min(cands, key=lambda k: E[k].max() - E[k][0])
    barrier = {s: (E[mep[s]].max() - E[mep[s]][0]) / N_FU * 1000 for s in mep}
    # Collective path in AlBN: the converged run whose saddle has every cation in the planar (layered-hexagonal) site
    col = ('AlBN', 'NEB')
    top = int(np.argmax(E[col])) + 1
    img = buy({'kind': 'neb_image', 'system': 'AlBN', 'run': 'NEB', 'image': top})
    off_col, _ = cation_offsets(img['cell_A'], img['species'], img['positions_A'])
    assert converged[col] and np.abs(off_col).max() < 0.2, off_col
    collective = (E[col].max() - E[col][0]) / N_FU * 1000
    # Born charges -> species averages for branch estimates
    zrows = buy({'kind': 'born_charges', 'system': 'AlBN'})['rows']
    zav = {sp: float(np.mean([r['Z33_e'] for r in zrows if r['atom'] == sp])) for sp in ('Al', 'B', 'N')}

    def berry_chain(system, folders):
        res = {f: buy({'kind': 'berry_phase', 'system': system, 'folder': '%02d_image' % f}) for f in folders}
        e_mep = E[mep[system]]
        for f, r in res.items():   # folder f is NEB image f+1 of the minimum-energy path (energies agree)
            assert abs(r['scf_total_energy_Ry'] * RY - e_mep[f]) < 5e-3, (system, f)
        vol = abs(np.linalg.det(np.array(res[folders[0]]['cell_A'])))
        p_acc = {folders[0]: 0.0}
        for a, b in zip(folders[:-1], folders[1:]):
            pa, q = berry_p(res[a])
            pb, _ = berry_p(res[b])
            dz = np.array(res[b]['positions_A'])[:, 2] - np.array(res[a]['positions_A'])[:, 2]
            guess = sum(zav[x] * d for x, d in zip(res[a]['species'], dz)) * E_CHARGE * 1e20 / vol
            d = pb - pa
            d += q * np.round((guess - d) / q)
            assert abs(d - guess) < 0.3 * q, (system, a, b, d, guess)
            p_acc[b] = p_acc[a] + d
        return res, p_acc

    # AlN: symmetric path; the midpoint image is the layered-hexagonal reference
    e_aln = E[mep['AlN']]
    mid = int(np.argmax(e_aln))
    res_n, acc_n = berry_chain('AlN', [0, mid // 2, mid])
    off_mid, _ = cation_offsets(res_n[mid]['cell_A'], res_n[mid]['species'], res_n[mid]['positions_A'])
    assert np.abs(off_mid).max() < 0.05
    ps_aln = -acc_n[mid]
    # AlBN: no layered-hexagonal image on the path; end states are inverted copies -> Ps = -dP(total)/2
    e_b = E[mep['AlBN']]
    saddle = int(np.argmax(e_b))
    minima = [k for k in range(1, len(e_b) - 1) if e_b[k] < e_b[k - 1] and e_b[k] < e_b[k + 1]]
    inter = min(minima, key=lambda k: e_b[k])
    last = len(e_b) - 1
    folders = sorted(set([0, 10, 20, saddle, inter, 50, 60, last]))
    res_b, acc_b = berry_chain('AlBN', folders)
    ps_albn = -acc_b[last] / 2
    saddle_dp = acc_b[saddle]
    off0, sp = cation_offsets(res_b[0]['cell_A'], res_b[0]['species'], res_b[0]['positions_A'])
    offi, _ = cation_offsets(res_b[inter]['cell_A'], res_b[inter]['species'], res_b[inter]['positions_A'])
    offs, _ = cation_offsets(res_b[saddle]['cell_A'], res_b[saddle]['species'], res_b[saddle]['positions_A'])
    inverted = int(np.sum(np.sign(offi) != np.sign(off0)))
    b_index = sp.index('B')
    inter_e = (e_b[inter] - e_b[0]) / N_FU * 1000
    mode = 'non-collective' if 0 < inverted < len(off0) and barrier['AlBN'] < collective else 'collective'
    submission = {'barrier_meV_per_fu': barrier, 'AlBN_collective_barrier_meV_per_fu': collective,
                  'spontaneous_polarization_C_m2': {'AlN': ps_aln, 'AlBN': ps_albn},
                  'AlBN_saddle_polarization_change_C_m2': saddle_dp,
                  'AlBN_intermediate': {'energy_meV_per_fu': inter_e, 'inverted_cations': inverted},
                  'AlBN_switching_mode': mode, 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    _, q_n = berry_p(res_n[0])
    _, q_b = berry_p(res_b[0])
    (APP/'report.md').write_text(f"""# Polarization switching in wurtzite AlN and (Al,B)N (reference solution)

**Runs.** neb.x runs were judged by their logs: AlN refined_NEB2 and AlBN NEB / refined_NEB2 report convergence
with all image errors below path_thr = 0.05 eV/A; AlBN refined_NEB stopped without convergence (killed job,
image errors up to {max(x['error_eV_per_A'] for x in neb[('AlBN', 'refined_NEB')]['final_table']):.2f} eV/A). No run
used a climbing image, so the barrier is the highest image; the paths are dense (25 and 75 images).
The supercells hold 32 atoms = 16 formula units.

**Barriers (per f.u.).** AlN {barrier['AlN']:.1f} meV (single symmetric hump; the saddle image has every cation
in its basal N plane: the layered-hexagonal state, collective switching). (Al,B)N, converged first-stage run:
{collective:.1f} meV, also collective (all cations planar at its saddle, max |offset| {np.abs(off_col).max():.2f} A):
only ~{(1 - collective / barrier['AlN']) * 100:.0f}% below AlN, which bounds what homogeneous effects (composition and the
~1% smaller lattice) do. The minimum-energy path (refined_NEB2) is {barrier['AlBN']:.1f} meV
(ratio to AlN {barrier['AlBN'] / barrier['AlN']:.2f}): switching is non-collective. Columns of cations invert in sequence;
the lowest intermediate minimum (image {inter + 1}, {inter_e:.1f} meV/f.u.) has {inverted} of 16 cations inverted, and at the
saddle (image {saddle + 1}) the B-containing column is passing through the planar site (B offset {offs[b_index]:+.2f} A).
The B column is not the first to invert. The archive has no non-collective AlN path and no AlN at the (Al,B)N
lattice, so neither B's specific role nor a pure strain effect is isolated by these data.

**Polarization.** pw.x prints P modulo the quantum eR/Omega of the supercell ({q_n:.4f} C/m^2 for AlN, {q_b:.4f} for
(Al,B)N). Branches were fixed between purchased images with Born-charge estimates (archived Z33, species
averages), using segments of <0.7 C/m^2 so the estimate error (~10-15%) is far below half a quantum; an
end-to-end estimate would pick the wrong branch. Berry folders NN_image are 0-based and match NEB image NN+1
(scf energies agree). Reference: the nonpolar layered-hexagonal structure (AlN: the centrosymmetric saddle
image, P = 0; (Al,B)N: midway between the inverted end states, Ps = -dP/2). Ps(AlN) = {ps_aln:.4f},
Ps(AlBN) = {ps_albn:.4f} C/m^2. A zinc-blende reference or raw end-state values (~0.1 C/m^2) would differ by more
than 1 C/m^2. At the (Al,B)N saddle the polarization has changed by {saddle_dp:.3f} C/m^2 (~{saddle_dp / (-2 * ps_albn) * 100:.0f}% of the reversal).

**Coercive field.** A homogeneous estimate E_c ~ Delta E/(Ps * Omega_fu) gives tens of MV/cm, far above measured
fields: real switching is nucleation and domain-wall limited, so only the trend (lower barrier, lower field in
(Al,B)N) is supported, not the magnitude.
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
