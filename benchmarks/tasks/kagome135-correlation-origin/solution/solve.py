"""Reference solution for kagome135-correlation-origin; runs as the Harbor oracle and in check_task.

Uses no literature values: electron counts from the archived Fermi levels and the staged models, model identity
from the .wout lattices, orbital names from the .wout projection table, all Z from lab evaluations.
"""
import json
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
MODELS = ['TiSb', 'VSb', 'CrSb', 'TiBi', 'CrBi', 'CrBi_9.2']
D_NAMES = {1: 'dz2', 2: 'dxz', 3: 'dyz', 4: 'dx2-y2', 5: 'dxy'}       # Wannier90 l = 2 real harmonics (mr)


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    result = json.loads(subprocess.run(cmd, capture_output=True, text=True).stdout)
    if 'error' in result:
        raise RuntimeError(result['error'])
    return result


def evaluate(request):
    r = lab('evaluate', request)
    while r['status'] == 'running':
        time.sleep(5)
        r = lab('job/%d' % r['id'])
    if r['status'] != 'complete':
        raise RuntimeError(str(r))
    return r['id'], r['result']


def wout(folder):
    return next((APP/'archive'/folder).glob('*.wout')).read_text()


def model(folder):
    d = np.load(APP/'models'/(folder + '.npz'))
    return d['R'].astype(float), d['degen'].astype(float), d['H']


def electron_count(folder, n=(12, 12, 8), kT=0.002):
    R, deg, H = model(folder)
    ef = float((APP/'archive'/folder/'Fermi_level.dat').read_text().split()[-2])
    g = np.stack(np.meshgrid(*[np.arange(m) / m for m in n], indexing='ij'), -1).reshape(-1, 3)
    total = 0.0
    for chunk in np.array_split(g, 8):
        hk = np.einsum('kr,rmn->kmn', np.exp(2j * np.pi * chunk @ R.T) / deg, H, optimize=True)
        e = np.linalg.eigvalsh(hk)
        total += (0.5 * (1 - np.tanh(0.5 * (e - ef) / kT))).sum()
    return 2 * total / len(g)


def sigma(folder):
    R, deg, H = model(folder)
    s2 = np.einsum('r,rab->a', 1 / deg, np.abs(H) ** 2)
    i0 = np.where((R == 0).all(1))[0][0]
    return np.sqrt(s2 - np.abs(np.diag(H[i0])) ** 2)


def main():
    # Pristine fillings: count the archive's own Fermi level (the charge-neutral DFT) in each 33-orbital model.
    raw = {f: electron_count(f) for f in MODELS}
    counts = {f: int(round(v)) for f, v in raw.items()}
    scf = (APP/'archive/DFT/scf.out').read_text()
    n_dft = float(re.search(r'number of electrons\s*=\s*(\S+)', scf).group(1))
    # Orbital names: first five projections (site 1) of the .wout table; l = 2 and the mr code.
    table = re.findall(r'\|\s+\S+\s+\S+\s+\S+\s+(\d)\s+(\d)\s+1\s', wout('CrSb'))
    names = [D_NAMES[int(mr)] for l, mr in table[:5]]
    assert all(l == '2' for l, mr in table[:15])
    # Which CsCr3Bi5 folder has the CsCr3Sb5 lattice?
    lat = {f: float(re.search(r'a_1\s+(\S+)', wout(f)).group(1)) for f in ('CrSb', 'TiBi', 'CrBi', 'CrBi_9.2')}
    sb_lattice = min(('CrBi', 'CrBi_9.2'), key=lambda f: abs(lat[f] - lat['CrSb']))
    tibi_lattice = 'CrBi_9.2' if sb_lattice == 'CrBi' else 'CrBi'
    evidence, results = [], {}

    def ssmf(f, n, U, ju=0.18):
        i, r = evaluate({'kind': 'ssmf', 'model': f, 'electrons': n, 'U_eV': U, 'J_over_U': ju})
        evidence.append(i)
        results[(f, n, U, ju)] = r
        return r

    def zs(r):
        return np.array([r['Z_by_wannier_function'][str(i)] for i in range(1, 6)])

    nT, nC = counts['TiSb'], counts['CrSb']
    cr = zs(ssmf('CrSb', nC, 6.0))
    ti = zs(ssmf('TiSb', nT, 6.0))
    ti_at_cr = zs(ssmf('TiSb', nC, 6.0))
    cr_at_ti = zs(ssmf('CrSb', nT, 6.0))
    cr_j0 = zs(ssmf('CrSb', nC, 6.0, 0.0))
    bi_sb = zs(ssmf(sb_lattice, counts[sb_lattice], 4.5))
    bi_tibi = zs(ssmf(tibi_lattice, counts[tibi_lattice], 4.5))
    # Context for the report (not graded): Sb at the same U, Bi at U = 6, the other pristine compounds,
    # the nominal-count variant, J scan.
    cr45 = zs(ssmf('CrSb', nC, 4.5))
    bi6 = ssmf(sb_lattice, counts[sb_lattice], 6.0)
    v = zs(ssmf('VSb', counts['VSb'], 6.0))
    tibi = zs(ssmf('TiBi', counts['TiBi'], 6.0))
    cr_nom = zs(ssmf('CrSb', nC - 1, 6.0))
    cr_j09 = zs(ssmf('CrSb', nC, 6.0, 0.09))
    L = lambda z: np.log(z.mean())
    total = L(cr) - L(ti)
    fill = 0.5 * ((L(ti_at_cr) - L(ti)) + (L(cr) - L(cr_at_ti)))
    imin = int(np.argmin(cr))
    sig_bi, sig_ti = sigma(sb_lattice)[imin], sigma('TiSb')[imin]
    submission = {
        'electron_counts': counts,
        'CsCr3Sb5_U6': {'Z_mean': float(cr.mean()), 'Z_min': float(cr.min())},
        'most_renormalized_orbital': names[imin],
        'filling_share': float(fill / total),
        'hund_ratio': float(cr_j0.mean() / cr.mean()),
        'CsCr3Bi5_U4p5': {'Sb_lattice_Z_mean': float(bi_sb.mean()), 'TiBi_lattice_Z_mean': float(bi_tibi.mean())},
        'sigma_eV': {'CsCr3Bi5_Sb_lattice': float(sig_bi), 'CsTi3Sb5': float(sig_ti)},
        'evidence_ids': evidence,
    }
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    fmt = lambda z: ', '.join('%s %.3f' % (names[i], z[i]) for i in range(5))
    (APP/'report.md').write_text(f"""# Correlations in 135 kagome metals (reference solution)

**Electron counts.** Counting the bands of each 33-orbital model below its archived Fermi level gives
{ {f: round(x, 2) for f, x in raw.items()} }, i.e. {counts}. The example DFT output has {n_dft:.0f} valence
electrons for CsV3Sb5, of which 46 bands (Cs 5s5p, V 3s3p, Sb 4d and 5s) are excluded from the Wannier window,
leaving 31: the Cs 6s electron is donated into the d-p manifold. Nominal counts that ignore it (27/30/33) are
one electron short; at that count CsCr3Sb5 gives Zbar = {cr_nom.mean():.3f} instead of {cr.mean():.3f}.

**Models.** The .wout lattices give a = {lat['CrBi']:.4f} A (CrBi) and {lat['CrBi_9.2']:.4f} A (CrBi_9.2), against
{lat['CrSb']:.4f} A for CrSb and {lat['TiBi']:.4f} A for TiBi: {sb_lattice} is CsCr3Bi5 at the CsCr3Sb5 lattice,
{tibi_lattice} at the CsTi3Bi5 lattice. The suffix '9.2' fits both c values (9.24 and 9.21 A) and does not
identify the lattice; the in-plane a does.

**Orbitals.** The .wout projection table lists, per site, l = 2 with mr = 5, 3, 1, 2, 4, i.e. {names} in the
site-local frame (z along c, x along the given axis, which points from the hexagon centre to the site); the
onsite Hamiltonian mixes only functions 3 and 5 (the two Ag orbitals dz2 and dx2-y2), confirming the order.
This is not Wannier90's default d order.

**Pristine CsCr3Sb5 (U = 6 eV, J = 0.18 U).** Z = {fmt(cr)}; Zbar = {cr.mean():.3f}, Zmin = {cr.min():.3f} in
{names[imin]}, the orbital with the smallest kinetic scale sigma (exact second moment of the PDOS,
sigma^2 = sum of |H_ab(R)|^2 over all hoppings out of the orbital). CsTi3Sb5: Z = {fmt(ti)}; CsV3Sb5:
Zbar = {v.mean():.3f}; CsTi3Bi5: Zbar = {tibi.mean():.3f}.

**Filling versus kinetic energy.** ln Zbar changes by {total:.3f} from CsTi3Sb5 (N = {nT}) to CsCr3Sb5
(N = {nC}). The Ti Hamiltonian at N = {nC} gives Zbar = {ti_at_cr.mean():.3f} (filling alone: almost nothing);
the Cr Hamiltonian at N = {nT} gives {cr_at_ti.mean():.3f} (Hamiltonian alone: modest). The two orders give very
different attributions; their average assigns {fill / total:.2f} of the change to the filling. The effect is a
synergy: filling matters only with the narrow Cr d bands, and the narrow bands matter mostly near half filling.
A cross-compound regression of Z on bandwidth and filling would confound the two.

**Hund versus Mott.** At U = 6 eV, J = 0 gives Zbar = {cr_j0.mean():.3f}, J = 0.09 U gives {cr_j09.mean():.3f},
J = 0.18 U gives {cr.mean():.3f}: the Hund coupling, not U, makes CsCr3Sb5 strongly correlated
(ratio {cr_j0.mean() / cr.mean():.2f}), with d occupations close to half filling per orbital. This is Hund's-metal
(and Hund-assisted orbital-selective Mott) physics; a U scan alone would not reveal it.

**Bi substitution.** At U = 4.5 eV: CsCr3Sb5 Zbar = {cr45.mean():.3f}; CsCr3Bi5 at the same lattice
{bi_sb.mean():.3f} (chemistry: Bi p states reduce the d kinetic scale, sigma({names[imin]}) = {sig_bi:.3f} eV vs
{sig_ti:.3f} eV in CsTi3Sb5); at the CsTi3Bi5 lattice {bi_tibi.mean():.3f} (lattice expansion adds a further,
larger reduction). These models have no spin-orbit coupling, so SOC plays no role in the effect. At U = 6 eV
CsCr3Bi5 (CsCr3Sb5 lattice) has metallic_solution = {bi6['metallic_solution']}: the metallic branch does not
survive, so it is the most correlated member of the family, on or beyond the verge of (orbital-selective) Mott
localization.

**Not supported.** Absolute Z values depend on the SSMF approximation (density-density interaction, FLL double
counting, mean-field slave spins) and on U and J, which are not computed here; the precise Mott boundary is
method dependent; no structural information on CsCr3Bi5 exists, so both lattices are assumptions; SOC, magnetism,
the CDW and dynamical effects are absent from these models. Comparison with ARPES or optics is not possible from
these data.
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
