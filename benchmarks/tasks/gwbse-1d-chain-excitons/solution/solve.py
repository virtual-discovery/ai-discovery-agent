"""Reference solution for gwbse-1d-chain-excitons; runs as the Harbor oracle and in check_task."""
import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np
from scipy.special import k0

APP = Path(os.environ.get('APP_DIR', '/app'))
CATALOG = json.loads((APP/'catalog.json').read_text())
HA, BOHR = 27.211386, 0.529177
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
        r = lab(f"job/{r['id']}")
    if r['status'] != 'complete':
        raise RuntimeError(r.get('error'))
    EVIDENCE.append(r['id'])
    return r['result']


def vbm(chain):
    return CATALOG['chains'][chain]['valence_band_max_index']


def qp_bands(result, chain, qp=True):
    """{k: (E_v, E_c)} from a qp purchase (QP energies, or DFT Eo when qp=False)."""
    out = {}
    v = vbm(chain)
    for k, rows in result['rows_by_k_index'].items():
        e = {b: eo + (de if qp else 0.0) for b, eo, de, sc in rows}
        out[int(k)] = (e[v], e[v + 1])
    return out


def direct_gap(result, chain):
    gaps = {k: c - v for k, (v, c) in qp_bands(result, chain).items()}
    k = min(gaps, key=gaps.get)
    return gaps[k], k


def min_direct_dft(bands):
    e = np.array(bands['energies_eV'])
    h = bands['highest_occupied_position']
    g = e[h + 1] - e[h]
    return float(g.min()), np.array(bands['k_frac']), g


def peaks(E, y, fraction):
    idx = [i for i in range(1, len(E) - 1) if y[i] > y[i - 1] and y[i] >= y[i + 1] and y[i] >= fraction * y.max()]
    return [float(E[i]) for i in idx]


def spectrum(result):
    d = np.array(result['data'])
    return d[:, 0], d[:, 1], d[:, 3]


def qp_database(result):
    line = next(l for l in result['header'] if 'KfnQPdb' in l)
    return line.split('"')[1]


def wire_binding(mu, alpha, r0, L=3000.0, N=1201):
    k = (np.arange(N) - N // 2) * 2 * np.pi / L
    q = np.abs(k[:, None] - k[None, :])
    q[q == 0] = 1e-12
    v = 2 * k0(q * r0)
    W = v / (1 + alpha * q ** 2 * v)
    h = 2 * np.pi / L
    W[np.arange(N), np.arange(N)] = 2 * (-np.log(h * r0 / 4) - 0.5772156649 + 1)
    return float(-np.linalg.eigvalsh(np.diag(k ** 2 / (2 * mu)) - W * h / (2 * np.pi))[0] * HA)


def main():
    # 1. SOC: archived DFT band structures with and without SOC
    bands = {(c, s): buy({'kind': 'dft_bands', 'chain': c, 'soc': s}) for c in ('Te3', 'Bi2Te3') for s in (True, False)}
    dsoc = {c: min_direct_dft(bands[(c, True)])[0] - min_direct_dft(bands[(c, False)])[0] for c in ('Te3', 'Bi2Te3')}
    # 2. QP data near the direct-gap minima (located from the cheap DFT bands)
    qp = {('Te3', 'scGW'): buy({'kind': 'qp', 'chain': 'Te3', 'level': 'scGW', 'k_indices': [20, 21, 22]}),
          ('Te3', 'G0W0'): buy({'kind': 'qp', 'chain': 'Te3', 'level': 'G0W0', 'k_indices': [21]}),
          ('S3', 'G0W0'): buy({'kind': 'qp', 'chain': 'S3', 'level': 'G0W0', 'k_indices': [22, 23]}),
          ('S3', 'scGW'): buy({'kind': 'qp', 'chain': 'S3', 'level': 'scGW', 'k_indices': [22, 23]}),
          ('As2S3', 'G0W0'): buy({'kind': 'qp', 'chain': 'As2S3', 'level': 'G0W0', 'k_indices': [1]}),
          ('As2S3', 'scGW'): buy({'kind': 'qp', 'chain': 'As2S3', 'level': 'scGW', 'k_indices': [1]}),
          ('Bi2Te3', 'G0W0'): buy({'kind': 'qp', 'chain': 'Bi2Te3', 'level': 'G0W0', 'k_indices': [1, 2, 3]})}
    gap = {key: direct_gap(r, key[0])[0] for key, r in qp.items()}
    # The GW starting point is the SOC band structure: compare Eo gaps with both DFT band structures.
    soc_check = {}
    for chain, key in (('Te3', ('Te3', 'scGW')), ('Bi2Te3', ('Bi2Te3', 'G0W0'))):
        eo = {k: c - v for k, (v, c) in qp_bands(qp[key], chain, qp=False).items()}
        for s in (True, False):
            _, kf, g = min_direct_dft(bands[(chain, s)])
            dev = max(abs(np.interp((k - 1) / 48, kf, g) - val) for k, val in eo.items())
            soc_check[f'{chain} soc={s}'] = round(float(dev), 4)
    assert soc_check['Te3 soc=True'] < 0.01 and soc_check['Te3 soc=False'] > 0.1
    # 3. BSE spectra; which QP energies each run used (header database and IP spectrum)
    runs = [('S3', 'BSE_G0W0'), ('S3', 'BSE_scGW'), ('Te3', 'BSE_G0W0'), ('Te3', 'BSE_scGW'),
            ('As2S3', 'BSE_G0W0'), ('As2S3', 'BSE_scGW'), ('Bi2Te3', 'BSE_G0W0')]
    bse = {key: buy({'kind': 'bse_spectrum', 'chain': key[0], 'run': key[1]}) for key in runs}
    db = {key: qp_database(r) for key, r in bse.items()}
    first_bright = {key: peaks(*spectrum(r)[:2], 0.10)[0] for key, r in bse.items()}
    ip_peak = {key: peaks(spectrum(r)[0], spectrum(r)[2], 0.01)[0] for key, r in bse.items()}
    same_energies = {}
    for chain in ('S3', 'Te3', 'As2S3'):
        a, b = spectrum(bse[(chain, 'BSE_G0W0')])[2], spectrum(bse[(chain, 'BSE_scGW')])[2]
        same_energies[chain] = bool(np.abs(a - b).max() < 1e-4 * b.max())
        assert same_energies[chain] == ('G0W0' not in db[(chain, 'BSE_G0W0')])
    optical = {c: first_bright[(c, 'BSE_scGW')] for c in ('S3', 'Te3', 'As2S3')}
    sc_shift = {c: gap[(c, 'scGW')] - gap[(c, 'G0W0')] for c in ('S3', 'Te3', 'As2S3')}
    dw = {}
    for chain in ('S3', 'Te3', 'As2S3'):
        raw = first_bright[(chain, 'BSE_scGW')] - first_bright[(chain, 'BSE_G0W0')]
        dw[chain] = raw if same_energies[chain] else raw - sc_shift[chain]
    # 4. Bi2Te3: its BSE_G0W0 run was computed on the scGW (G3W3) QP energies, so its IP spectrum gives the
    #    scGW transitions; the transferred self-consistency increment of the other chains does not apply.
    assert 'G3W3' in db[('Bi2Te3', 'BSE_G0W0')]
    bi_gap = ip_peak[('Bi2Te3', 'BSE_G0W0')]
    bi_sc_shift = bi_gap - gap[('Bi2Te3', 'G0W0')]
    gaps = {'Te3': {'with_soc': gap[('Te3', 'scGW')], 'without_soc': gap[('Te3', 'scGW')] - dsoc['Te3']},
            'Bi2Te3': {'with_soc': bi_gap, 'without_soc': bi_gap - dsoc['Bi2Te3']}}
    # Screened-interaction update for Bi2Te3: scale its W0 binding energy by the released chains' ratio,
    # taken proportional to the relative self-consistent gap opening that drives the change of W.
    eb_w0 = bi_gap - first_bright[('Bi2Te3', 'BSE_G0W0')]
    slopes = []
    for chain in ('S3', 'Te3', 'As2S3'):
        eb_sc = gap[(chain, 'scGW')] - optical[chain]
        slopes.append((eb_sc / (eb_sc + dw[chain]) - 1) / (sc_shift[chain] / gap[(chain, 'G0W0')]))
    eb_bi = eb_w0 * (1 + float(np.mean(slopes)) * bi_sc_shift / gap[('Bi2Te3', 'G0W0')])
    # 5. A 1D screened-wire model with the DFPT polarizability, for comparison only
    scr = buy({'kind': 'wire_screening', 'chain': 'Bi2Te3'})
    cell = np.array(scr['dfpt_cell_A'])
    alpha = (scr['epsilon_infinity_supercell'][2][2] - 1) * abs(np.linalg.det(cell[:2, :2])) / BOHR ** 2 / (4 * np.pi)
    c_bohr = scr['chain_lattice_constant_A'] / BOHR
    r0 = float(np.sqrt(scr['environ_pv_energy_Ry'] / scr['environ_pressure_Ry_per_bohr3'] / c_bohr / np.pi))
    e = qp_bands(qp[('Bi2Te3', 'G0W0')], 'Bi2Te3')
    dk = 2 * np.pi / (48 * c_bohr)
    curv = [(e[1][j] - 2 * e[2][j] + e[3][j]) / HA / dk ** 2 for j in (0, 1)]
    mu = 1 / (abs(curv[1]) + abs(curv[0]))
    eb_wire = wire_binding(mu, alpha, r0)

    submission = {'scgw_direct_gap_eV': gaps, 'scgw_optical_gap_eV': optical, 'screening_update_shift_eV': dw,
                  'bi2te3_binding_energy_eV': eb_bi, 'evidence_ids': EVIDENCE}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# Exfoliable 1D chains from the archived GW-BSE data (reference solution)

**SOC.** The Eo energies of the archived QP tables reproduce the SOC band structures: on the purchased k-points
the direct DFT gap deviates from the SOC bands by at most {max(soc_check['Te3 soc=True'], soc_check['Bi2Te3 soc=True']):.3f} eV
and from the SOC-free bands by up to {max(soc_check['Te3 soc=False'], soc_check['Bi2Te3 soc=False']):.2f} eV. The paired spinor bands are split: the GW and BSE runs already include
spin-orbit coupling, so no SOC correction is added to them. The header line "Perturbative SOC: no" only means
that SOC was not added perturbatively. SOC lowers the minimum direct DFT gap by {-dsoc['Te3']:.3f} eV (Te3) and
{-dsoc['Bi2Te3']:.3f} eV (Bi2Te3); the SOC-free gaps add these back, assuming SOC-independent QP corrections.

**Run labels.** The KfnQPdb lines show that the BSE_G0W0 runs of Te3, As2S3 and Bi2Te3 used the scGW QP
energies ({'; '.join(f'{c} {r}: {v}' for (c, r), v in db.items())}); their independent-particle spectra coincide with those of the BSE_scGW runs, so the two runs
differ only in the screened interaction of the kernel. Only S3's BSE_G0W0 run used G0W0 energies. The screening
update shifts (scGW-level W minus W0 at scGW energies) are {', '.join(f'{c} {v:+.3f}' for c, v in dw.items())} eV;
for S3 the QP shift of the minimum direct gap ({sc_shift['S3']:.3f} eV) was removed first.

**Self-consistency.** Eigenvalue self-consistency opens the direct gap by
{', '.join(f'{c} {v:.3f}' for c, v in sc_shift.items())} eV, but only by {bi_sc_shift:.3f} eV for Bi2Te3,
whose scGW gap ({bi_gap:.3f} eV) is read from the IP onset of its BSE_G0W0 run (scGW energies). Transferring the
mean increment of the other chains would put it near {gap[('Bi2Te3', 'G0W0')] + np.mean(list(sc_shift.values())):.2f} eV.

**Excitons.** First bright excitons (10% rule) at the scGW level: {', '.join(f'{c} {v:.3f}' for c, v in optical.items())}
eV. Binding energies relative to the minimum direct QP gap are 2.1-2.3 eV for S3 and Te3: the bare Coulomb
interaction of a thin wire is weakly screened (the supercell epsilon_zz only reflects the vacuum fraction;
the relevant quantity is the polarizability per length, {alpha:.0f} bohr^2 for Bi2Te3). QP gap opening and
binding largely cancel, so optical gaps lie only 0.4-1.5 eV above the DFT gaps. As2S3's band-edge transition at
Gamma is nearly dark along the chain; its first bright exciton belongs to higher transitions, and weak features
below it carry 1-2% of the peak height. Te3 and S3 band-edge transitions are also weak. Updating W to the scGW
level always strengthens binding.

**Bi2Te3 transfer.** Its BSE_G0W0 run (scGW energies, W0) gives E_b(W0) = {eb_w0:.3f} eV; scaling by the W-update
ratio of the released chains, taken proportional to the relative self-consistent gap opening, gives E_b =
{eb_bi:.3f} eV. A 1D screened-wire model with the DFPT polarizability (r0 = {r0:.1f} bohr from the quantum volume,
mu = {mu:.2f}) gives {eb_wire:.2f} eV and overestimates the binding of all chains, so it does not transfer.

**Not supported.** k-point, vacuum and band convergence (one grid, one cutoff radius; Bi2Te3 used 700 bands
versus 1450-1600 for the others); G0W0+BSE optical gaps for Te3, As2S3 and Bi2Te3 (not archived); the energies of
dark excitons (Haydock spectra only show bright ones); the exact Bi2Te3 W-update correction.
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
