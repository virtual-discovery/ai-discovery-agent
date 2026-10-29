"""Reference solution for heis-configdesign-fe; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


# bcc shells of the spiral cell (a = 2.87 A); shell radii scale with the lattice constant.
A = np.array(CATALOG['spiral_cell']['primitive_cell_A'])
A0 = CATALOG['spiral_cell']['lattice_constant_A']
N = np.array([(i, j, k) for i in range(-7, 8) for j in range(-7, 8) for k in range(-7, 8) if (i, j, k) != (0, 0, 0)])
D = np.round(np.linalg.norm(N @ A, axis=1), 4)
RADII = sorted(set(D))[:10]
B = 2 * np.pi * np.linalg.inv(A).T


def spiral_row(q, cone_deg, nshell):
    s2 = np.sin(np.radians(cone_deg)) ** 2
    return np.array([(1 - np.cos(2 * np.pi * N[D == r] @ np.array(q))).sum() * s2 for r in RADII[:nshell]])


def random_row(cell, frac, spins, nshell):
    cell, frac, e = np.array(cell), np.array(frac), np.array(spins, float)
    e = e / np.linalg.norm(e, axis=1, keepdims=True)
    pos = frac @ cell
    near = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) @ cell
    d0 = np.linalg.norm(pos[None] + near[:, None] - pos[0], axis=2)
    a = d0[d0 > 1e-6].min() * 2 / np.sqrt(3)          # this supercell's own bcc lattice constant
    radii = [r * a / A0 for r in RADII[:nshell]]
    span = int(np.ceil(radii[-1] / np.linalg.norm(cell, axis=1).min())) + 1
    shifts = np.array([(i, j, k) for i in range(-span, span + 1) for j in range(-span, span + 1)
                       for k in range(-span, span + 1)]) @ cell
    row = np.zeros(nshell)
    for i in range(len(pos)):
        dist = np.linalg.norm(pos[None] + shifts[:, None] - pos[i], axis=2)
        for s, r in enumerate(radii):
            row[s] += ((np.abs(dist - r) < 2e-3) * (1 - e @ e[i])[None, :]).sum()
    return row / len(pos)


def fit(X, y):
    coef = np.linalg.lstsq(np.c_[np.ones(len(y)), X], y, rcond=None)[0]
    return coef[0], coef[1:]


def main():
    evidence, spiral_rows, random_rows = [], [], []
    # Spirals: several non-[001] directions (single directions alias shells), plus small q for the stiffness.
    spread = [0, 1, 3, 5, 7, 10, 13, 16, 19]
    plan = {'[111]': spread, '[113]': spread, '[121]': spread, '[001]': [0, 1]}
    keys = {(s['direction'], s['q_index']): s['key'] for s in CATALOG['spirals_available']}
    for direction, indices in plan.items():
        for k in indices:
            r = lab('evaluate', {'kind': 'spiral', 'key': keys[(direction, k)]})
            evidence.append(r['id'])
            spiral_rows.append((direction, r['result']))
    # Random states in the larger supercells only (2x2x1 cells alias the first shells).
    labels = [s for s in CATALOG['random_states_available'] if s['supercell'] == '236']
    labels += [s for s in CATALOG['random_states_available'] if s['supercell'] == '341'][:10]
    for s in labels:
        r = lab('evaluate', {'kind': 'random_state', 'label': s['label']})
        evidence.append(r['id'])
        random_rows.append((s, r['result']))
    # Low-angle model (cone factor sin^2 theta); the shell count is chosen by leave-one-direction-out
    # cross-validation rather than by the in-sample error.
    def spiral_fit(rows, nshell):
        X = np.array([spiral_row(res['q_reciprocal_lattice'], res['cone_angle_deg'], nshell) for d, res in rows])
        y = np.array([res['energy_meV_per_atom'] for d, res in rows])
        return fit(np.r_[np.zeros((1, nshell)), X], np.r_[0.0, y])
    cv = {}
    for nshell in range(4, 11):
        errors = []
        for held in ('[111]', '[113]', '[121]'):
            a, Jh = spiral_fit([row for row in spiral_rows if row[0] != held], nshell)
            test = [res for d, res in spiral_rows if d == held]
            pred = [a + spiral_row(res['q_reciprocal_lattice'], res['cone_angle_deg'], nshell) @ Jh for res in test]
            errors.append(np.sqrt(np.mean([(p - res['energy_meV_per_atom']) ** 2 for p, res in zip(pred, test)])))
        cv[nshell] = float(np.mean(errors))
    nsp = min(cv, key=cv.get)
    c0, J = spiral_fit(spiral_rows, nsp)
    spiral_pred = {p['id']: float(c0 + spiral_row(p['q_reciprocal_lattice'], p['cone_angle_target_deg'], nsp) @ J)
                   for p in CATALOG['spirals_to_predict']}
    # Stiffness from the small-q limit: E / sin^2(theta) = D q^2, q in 1/A including 2 pi.
    vals = []
    for direction in plan:
        pts = sorted([res for d, res in spiral_rows if d == direction],
                     key=lambda res: np.linalg.norm(np.array(res['q_reciprocal_lattice']) @ B))[:2]
        q2 = np.array([np.linalg.norm(np.array(res['q_reciprocal_lattice']) @ B) ** 2 for res in pts])
        e = np.array([res['energy_meV_per_atom'] / np.sin(np.radians(res['cone_angle_deg'])) ** 2 for res in pts])
        vals.append((e * q2).sum() / (q2 ** 2).sum())
    stiffness = float(np.mean(vals))
    # Large-angle states: a separate effective model fitted to random states (the low-angle model does not transfer).
    nr = 6
    Xr = np.array([random_row(s['cell_A'], s['fractional_positions'], res['moments_muB'], nr) for s, res in random_rows])
    yr = np.array([res['energy_meV_per_atom'] for s, res in random_rows])
    r0, Jr = fit(Xr, yr)
    random_pred = {p['id']: float(r0 + random_row(p['cell_A'], p['fractional_positions'], p['spin_directions'], nr) @ Jr)
                   for p in CATALOG['random_states_to_predict']}
    low_angle_on_random = np.mean([random_row(p['cell_A'], p['fractional_positions'], p['spin_directions'], nsp) @ J
                                   for p in CATALOG['random_states_to_predict']])
    submission = {'spiral_predictions_meV_per_atom': spiral_pred, 'random_predictions_meV_per_atom': random_pred,
                  'J1_meV': float(J[0]), 'spin_wave_stiffness_meV_A2': stiffness, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# bcc Fe exchange model (reference solution)

Convention: E/atom - E_FM = (1/N) sum_i sum_j J_ij (1 - e_i.e_j), unit vectors, ordered pairs. The released
spirals are conical (cone angle ~30 deg), so each spiral energy carries a sin^2(theta) factor.
Low-angle model: {nsp} shells (chosen by leave-one-direction-out validation, errors {cv}) fitted to [111], [113]
and [121] spirals plus small-q [001] points; [001] alone aliases shells (J1 = J2) and single directions are rank
deficient. J1 = {J[0]:.2f} meV.
Stiffness from the two smallest q per direction: D = {stiffness:.1f} meV A^2 (flat spiral, q in 1/A with 2 pi).
Random states were computed at a = 2.83 A (spirals: 2.87 A); shells were assigned in each cell's own lattice.
The low-angle model predicts the hidden 3x3x4 random states at {low_angle_on_random:.0f} meV/atom on average, far
below the large-angle data: the Heisenberg picture fitted at small angles does not transfer to large-angle
disorder (non-Heisenberg/itinerant effects, with the volume difference as a confound). A separate {nr}-shell
effective model fitted to 2x3x6 and 3x4x1 states predicts them instead. Long-range tails, and anything that
depends on them (Curie temperature, truncated stiffness sums), are poorly constrained by these data.
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
