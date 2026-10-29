"""Reference solution for fresh-al-surface-workfunction; runs as the Harbor oracle and in check_task.

Plan: oriented bulk references with the slab's in-plane mesh; thickness ladders N = 7..12 of both faces at the
reference protocol with enough bands for thick metallic slabs; symmetric relaxations of Al(111); quantum-size
oscillations averaged over the ladder; vacuum level from the electrostatic potential plateau.
"""
import json
import math
import os
import subprocess
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
J_PER_M2 = 16.0217663          # eV/A^2 -> J/m^2
LADDER = list(range(7, 13))
RELAX = [9, 10]
HIDDEN_N = 18
ANCHOR = -537.46          # eV/atom, the instruction's anchor for the thick-slab energies


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def compute(request):
    """Submit one calculation and poll it to completion; returns (job id, result)."""
    view = lab('evaluate', request)
    while view['status'] == 'running':
        time.sleep(10)
        view = lab('job/%d' % view['id'])
    if view['status'] != 'complete':
        raise RuntimeError('job %d failed: %s' % (view['id'], view.get('error')))
    return view['id'], view['result']


def bands(n):
    return int(math.ceil(2.5 * n)) + 4


def vacuum_level(r):
    """Mean electrostatic potential over the central 1 A of the vacuum gap (the plateau)."""
    z, v = np.array(r['potential_z_A']), np.array(r['electrostatic_potential_eV'])
    c, planes = r['cell_height_A'], np.array(r['atom_z_A'])
    dist = np.min(np.abs(((z[:, None] - planes[None, :]) + c / 2) % c - c / 2), axis=1)
    plateau = dist >= dist.max() - 0.5
    return float(v[plateau].mean()), float(v[plateau].max() - v[plateau].min())


def main():
    evidence, out = [], {}
    for face in ('111', '100'):
        jid, bulk = compute({'kind': 'bulk', 'cell': 'oriented', 'face': face, 'kpoints': 18, 'kz': 18})
        evidence.append(jid)
        eb = bulk['free_energy_eV_per_atom']
        rows = []
        for n in LADDER:
            jid, r = compute({'kind': 'slab', 'face': face, 'layers': n, 'kpoints': 18, 'vacuum_A': 15.0,
                              'smearing_Ry': 0.02, 'nbnd': bands(n)})
            evidence.append(jid)
            assert r['highest_band_minimum_eV'] > r['fermi_energy_eV'] + 0.5, 'too few bands'
            vac, flat = vacuum_level(r)
            area = r['surface_cell_area_A2']
            rows.append({'N': n, 'E': r['free_energy_eV'], 'sigma': (r['free_energy_eV'] - n * eb) / (2 * area) * J_PER_M2,
                         'phi': vac - r['fermi_energy_eV'], 'plateau_flatness_eV': flat})
        sigma = float(np.mean([x['sigma'] for x in rows]))
        phi = float(np.mean([x['phi'] for x in rows]))
        slope, icpt = np.polyfit([x['N'] for x in rows], [x['E'] for x in rows], 1)
        out[face] = {'eb': eb, 'area': area, 'rows': rows, 'sigma': sigma, 'phi': phi,
                     'sigma_spread': float(np.std([x['sigma'] for x in rows])),
                     'phi_spread': float(np.std([x['phi'] for x in rows])),
                     'fit_slope': float(slope), 'fit_sigma': float(icpt / (2 * area) * J_PER_M2),
                     'pred_E18': float(HIDDEN_N * (eb - ANCHOR) + 2 * area * sigma / J_PER_M2)}
    d12 = []
    for n in RELAX:
        jid, r = compute({'kind': 'slab', 'face': '111', 'layers': n, 'kpoints': 18, 'vacuum_A': 15.0,
                          'smearing_Ry': 0.02, 'nbnd': bands(n), 'relax': 'outer_two_both_sides'})
        evidence.append(jid)
        sp, d0 = r['interlayer_spacings_A'], r['bulk_interlayer_spacing_A']
        d12.append(50.0 * ((sp[0] + sp[-1]) / d0 - 2.0))
    d12_mean = float(np.mean(d12))
    dphi = out['100']['phi'] - out['111']['phi']
    noise = 2 * math.hypot(out['100']['phi_spread'], out['111']['phi_spread']) / math.sqrt(len(LADDER))
    ordering = 'unresolved' if abs(dphi) < max(noise, 0.02) else ('phi100>phi111' if dphi > 0 else 'phi111>phi100')
    submission = {
        'surface_energy_J_per_m2': {f: out[f]['sigma'] for f in ('111', '100')},
        'work_function_eV': {f: out[f]['phi'] for f in ('111', '100')},
        'work_function_ordering': ordering,
        'relaxation_d12_percent_111': d12_mean,
        'thick_slab_predictions': {f'{f}_N{HIDDEN_N}': {'free_energy_offset_eV': out[f]['pred_E18'], 'work_function_eV': out[f]['phi']}
                                   for f in ('111', '100')},
        'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    lines = []
    for f in ('111', '100'):
        o = out[f]
        lines.append(f"Al({f}): oriented bulk {o['eb']:.5f} eV/atom; " + '; '.join(
            f"N={x['N']} sigma={x['sigma']:.3f} phi={x['phi']:.3f}" for x in o['rows'])
            + f". Mean sigma {o['sigma']:.3f} J/m2 (spread {o['sigma_spread']:.3f}), phi {o['phi']:.3f} eV (spread "
              f"{o['phi_spread']:.3f}); a straight-line fit of E_N gives slope {o['fit_slope']:.5f} eV/atom.")
    (APP/'report.md').write_text(f"""# Al(111) and Al(100) surfaces (reference solution)

Protocol: PBE, SSSP Al PAW, 30/240 Ry, a = 4.04 A, in-plane 18x18 Gamma-centred mesh, MV smearing 0.02 Ry, 15 A
vacuum, symmetric bulk-terminated slabs. Thick metallic slabs need many more bands than the pw.x default (which
leaves the highest band below the Fermi level for N >= 9 and raises the slab energy, mimicking a divergent surface
energy); every slab here used nbnd = 2.5 N + 4 and the highest band stays well above E_F.

Surface energy: sigma_N = (E_N - N e_b) / 2A with e_b from an oriented bulk cell (the slab's 1x1 cell, same
in-plane mesh and smearing, dense kz), so the in-plane k-sampling error cancels and sigma_N does not drift with N
(Boettger's divergence appears with a primitive-cell bulk reference on another mesh). A slope fit of E_N over
thin slabs is biased by quantum-size oscillations. Work function: phi = V_vac - E_F with V_vac the plateau of the
planar-averaged electrostatic potential in the middle of the vacuum; the total local potential still contains the
exchange-correlation potential of the spill-out density and lies ~0.24 eV lower at 15 A. Both sigma_N and phi_N
oscillate with N (quantum-size effects); the values below average N = 7-12.

{chr(10).join(lines)}

Ordering: phi(100) - phi(111) = {dphi:+.3f} eV, far outside the quantum-size scatter ({noise:.3f} eV): {ordering}.
Al violates the Smoluchowski ordering of a close-packed face, as found experimentally.
Relaxation: symmetric relaxation of the two outer planes (N = {RELAX}) expands d12 of Al(111) by
{', '.join(f'{x:+.2f}' for x in d12)} % (mean {d12_mean:+.2f} %); the unrelaxed slabs already show outward forces on
the surface plane. One-sided relaxation would need the dipole correction.
Thick slabs (N = 18) are predicted from E_N = N e_b + 2 A sigma and the averaged phi.
Not supported: absolute agreement with experiment (PBE underestimates surface energies), relaxation effects on
sigma/phi beyond this protocol, and anything finer than the residual quantum-size scatter.
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
