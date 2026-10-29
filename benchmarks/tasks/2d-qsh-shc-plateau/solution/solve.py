"""Reference solution for 2d-qsh-shc-plateau; runs as the Harbor oracle and in check_task."""
import json
import math
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
E2H = 3.874045864e-5            # e^2/h in S
CANTED = ('Bi2_P-3m1', 'Hg4O2_P1')


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(3)
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(f'job {r["id"]} failed: {r.get("error")}')
    return r['id'], r['result']


def archive_facts(name):
    d = APP/'archive'/name
    out, win, inp = (d/'aiida.out.excerpt').read_text(), (d/'aiida.win').read_text(), (d/'aiida.in').read_text()
    nelec = int(round(float(re.search(r'number of electrons\s*=\s*(\S+)', out).group(1))))
    excl = re.search(r'exclude_bands\s*=\s*([^\n]+)', win).group(1)
    nexcl = sum(int(b) - int(a) + 1 for a, b in re.findall(r'(\d+)\s*-\s*(\d+)', excl))
    cell = np.array([[float(x) for x in l.split()] for l in inp.split('CELL_PARAMETERS')[1].split('\n')[1:4]])
    return {'nocc': nelec - nexcl, 'c_A': abs(cell[2, 2]),
            'scf_fermi': float(re.search(r'the Fermi energy is\s+(\S+)', out).group(1)),
            'win_vbm': float(re.search(r'# VBM = (\S+)', win).group(1)),
            'win_cbm': float(re.search(r'# CBM = (\S+)', win).group(1))}


def z2_from_wcc(w):
    """Parity of the crossings of the largest-gap midpoint by the hybrid Wannier centres (k from 0 to 1/2)."""
    def ref(a):
        a = np.sort(a)
        g = np.diff(np.r_[a, a[0] + 1])
        j = int(np.argmax(g))
        return (a[j] + g[j] / 2) % 1
    total = 0
    for i in range(len(w) - 1):
        z0, z1 = ref(w[i]), ref(w[i + 1])
        b = np.sort(w[i + 1])
        lo, hi = min(z0, z1), max(z0, z1)
        total += int(np.sum((b > lo) & (b < hi))) if hi - lo < 0.5 else int(np.sum((b > hi) | (b < lo)))
    return total % 2


def main():
    evidence, notes = [], {}
    z2, plateau, direction = {}, {}, {}
    for name in CATALOG['materials']:
        f = archive_facts(name)
        # 1. Band edges of the Wannier model on a dense mesh: the scf Fermi level (cold smearing) and the .win VBM/CBM
        #    (coarse scf grid) are not reliable positions for the chemical potential.
        i, be = evaluate({'kind': 'band_edges', 'material': name, 'band': f['nocc'], 'mesh': 300})
        evidence.append(i)
        vbm, cbm = be['max_of_band_n']['energy_eV'], be['min_of_band_n_plus_1']['energy_eV']
        gap, mid = cbm - vbm, 0.5 * (vbm + cbm)
        # 2. Z2 from the Wilson-loop flow of the occupied bands; small-gap models need a fine k sampling.
        nk, nl = (201, 119) if gap < 0.1 else (51, 110)
        i, w = evaluate({'kind': 'wcc', 'material': name, 'bands': f['nocc'], 'loop_axis': 1,
                         'k': [round(x, 6) for x in np.linspace(0, 0.5, nk)], 'n_loop': nl})
        evidence.append(i)
        z2[name] = z2_from_wcc(np.array(w['wcc']))
        # 3. Spin Hall plateau at mid-gap, clean limit, Wannier-centre position operator (the lattice convention
        #    breaks the crystal symmetry of the response), with adaptive refinement for small gaps.
        spins = ['x', 'y', 'z'] if name in CANTED else ['z']
        req = {'kind': 'shc', 'material': name, 'spin': spins, 'eta_eV': 0.0005, 'position': 'wannier_centres',
               'mu_eV': [round(mid, 6)]}
        if gap > 0.2:
            req['mesh'] = [60, 60]
        else:
            req.update(mesh=[120, 120], refine={'factor': 8, 'threshold_A2': 10.0, 'max_points': 400})
        i, s = evaluate(req)
        evidence.append(i)
        conv = f['c_A'] * 1e-8 / E2H          # (hbar/e) S/cm -> (hbar/e)(e^2/h) per sheet, with this cell's height
        vec = {sp: s['sigma_xy_hbar_over_e_S_per_cm'][sp][0] * conv for sp in spins}
        plateau[name] = vec['z']
        if name in CANTED:
            x, y, z = vec['x'], vec['y'], vec['z']
            direction[name] = {'theta_deg': math.degrees(math.atan2(math.hypot(x, y), z)),
                               'phi_deg': math.degrees(math.atan2(y, x))}
        notes[name] = {'gap_eV': gap, 'mid_eV': mid, 'scf_fermi_minus_vbm': f['scf_fermi'] - vbm,
                       'win_mid_minus_vbm': 0.5 * (f['win_vbm'] + f['win_cbm']) - vbm, 'sheet': vec,
                       'refined': s['refined_points']}
        # 4. Convergence check for the smallest gap (k-mesh aliasing of the curvature hot spots).
        if gap < 0.03:
            req2 = dict(req, mesh=[180, 180])
            i, s2 = evaluate(req2)
            evidence.append(i)
            notes[name]['check_180'] = s2['sigma_xy_hbar_over_e_S_per_cm']['z'][0] * conv
    submission = {'z2': z2, 'plateau_sigma_z': plateau, 'spin_direction': direction, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    lines = [f"- {n}: Z2 = {z2[n]}, gap {notes[n]['gap_eV'] * 1e3:.1f} meV, mu = {notes[n]['mid_eV']:.4f} eV, "
             f"sigma^z = {plateau[n]:.3f} (hbar/e)(e^2/h); scf E_F - VBM = {notes[n]['scf_fermi_minus_vbm'] * 1e3:.0f} meV; "
             f".win mid-gap - VBM = {notes[n]['win_mid_minus_vbm'] * 1e3:.0f} meV" for n in plateau]
    (APP/'report.md').write_text(f"""# Spin Hall plateaus of eight archived 2D insulators (reference solution)

Occupied Wannier bands = electrons of the scf run minus `exclude_bands`. For each archived Wannier model the band
edges were located on a 300 x 300 mesh. The scf Fermi level (cold smearing) lies inside the valence band of the
small-gap models and the `.win` fermi_energy is the coarse-grid VBM; for Bi2In2 the `.win` VBM/CBM come from the scf
grid and the Wannier VBM at a generic k lies above the `.win` mid-gap. The chemical potential was therefore placed at
the mid-gap of the Wannier model.

Z2 comes from the Wilson-loop flow of the hybrid Wannier centres over half the Brillouin zone (largest-gap reference
line); small-gap models need a fine k sampling (coarse samplings miss the fast winding near the gap minimum).

Spin Hall conductivities were computed at T = 0 with eta = 0.5 meV (clean limit), with the velocity operator built
from the Wannier centres: placing all Wannier functions on the lattice points breaks the crystal symmetry of the
response (Bi2 then shows a sigma^x_xy forbidden by its mirror). Small-gap models were refined adaptively. The service
returns Wannier90 units with the vacuum-padded cell; values were converted per sheet with each cell height c,
sigma_2D = sigma_W90 * c / (e^2/h).

{chr(10).join(lines)}

CuLi2As (D3h, horizontal mirror) conserves S_z and is quantized. The other Z2-odd models do not conserve S_z, so the
plateau is not quantized and the spin of the spin current tilts away from z where symmetry allows in-plane components
(Bi2: sigma^y_xy allowed and sigma^x_xy forbidden by the mirror x -> -x; Hg4O2 (P1): all components allowed). Bi2In2 is
Z2-trivial yet has a large in-gap SHC: the SHC magnitude does not determine the topology. Spin directions:
{json.dumps(direction)}.

Limitations: the spin operator is the projected spinor operator (no .spn in the archive) and the position operator
is diagonal at the Wannier centres; the archived Top_materials SHC files are for relaxed structures with other lattice
constants and do not describe these models; DFT-PBE band gaps, and the phonon stability of several candidates is not
established by these data.
""")
    lab('submit', submission)
    for _ in range(120):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(5)


if __name__ == '__main__':
    main()
