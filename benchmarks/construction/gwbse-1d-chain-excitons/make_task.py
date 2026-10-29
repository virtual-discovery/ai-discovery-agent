"""Generate the staged, lab, private and ablation data for gwbse-1d-chain-excitons.

    python3 benchmarks/construction/gwbse-1d-chain-excitons/fetch.py
    cd benchmarks/construction/gwbse-1d-chain-excitons && python3 -B make_task.py

Reads the verified archive members from ~/.cache/ai-theorist/gwbse-1d-chain-excitons and writes
catalog.json, lab/data.json, private/hidden.json, private/hidden_values.json,
solution/ablation_values.json and the budget in both public.json files. Prints the calibration table.
"""
import json
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/gwbse-1d-chain-excitons'
BUDGET = 45
COSTS = {'qp_per_k_point': 1, 'bse_spectrum': 3, 'dft_bands': 1, 'wire_screening': 1}
WITHHELD = {'qp': [('Bi2Te3', 'scGW')], 'bse': [('Bi2Te3', 'BSE_scGW')]}
THRESHOLDS = {'gap_abs_eV': 0.05, 'optical_gap_abs_eV': 0.05, 'screening_shift_abs_eV': 0.10,
              'bi2te3_binding_window_eV': [0.68, 1.15]}
FORMULA_ATOMS = {'S3': 3, 'Te3': 3, 'As2S3': 5, 'Bi2Te3': 5}
SOC_WINDOW = (12, 12)     # occupied, empty bands served from the SOC band files
NOSOC_WINDOW = (6, 6)
SOURCE = 'retrieved author-archived result (no new calculation)'


def r(x, n=6):
    return float(round(float(x), n))


def qp_payload(chain, level):
    qp = M.read_qp(chain, level)
    rows = qp['rows']
    by_k = {}
    for k, b, eo, de, sc in rows:
        by_k.setdefault(str(int(k)), []).append([int(b), r(eo), r(de), r(sc)])
    return {'file': qp['file'], 'header': qp['header'], 'row_format': ['band', 'Eo [eV]', 'E-Eo [eV]', 'Sc|Eo [eV]'],
            'rows_by_k_index': by_k}


def bse_payload(chain, run):
    b = M.read_bse(chain, run)
    return {'file': b['file'], 'header': b['header'], 'columns': b['columns'],
            'data': [[float('%.6g' % x) for x in row] for row in b['data']]}


def bands_payload(chain, soc):
    b = M.read_bands(chain, soc)
    v = M.VBM[chain] if soc else M.VBM_NOSOC[chain]
    occ, emp = SOC_WINDOW if soc else NOSOC_WINDOW
    window = b['energies'][v - occ:v + emp]
    return {'file': b['file'], 'spin_orbit_coupling': soc, 'k_frac': [r(k, 5) for k in b['k_frac']],
            'energies_eV': [[r(e, 4) for e in band] for band in window],
            'highest_occupied_position': occ - 1,
            'note': (f'{occ} highest occupied and {emp} lowest empty bands of the archived file, ordered by energy; '
                     'entry highest_occupied_position is the highest occupied band. Absolute PBE eigenvalues as '
                     'archived (no Fermi level is archived). k_frac is k along Gamma-Z in units of 2*pi/c.')}


def screening_payload(chain):
    s = M.read_screening(chain)
    return {**s, 'chain_lattice_constant_A': M.CHAIN_C_A[chain],
            'note': ('epsilon_infinity_supercell: DFPT electronic dielectric tensor of the periodic supercell '
                     '(chain along z, vacuum in x and y), from IR/<chain>/ph.epsil.out. The Environ record is the '
                     'quantum-volume run of volume/<chain>: the PV term of the energy divided by the applied '
                     'pressure is the quantum volume of one cell.')}


def main():
    data = {'costs': COSTS, 'qp': {}, 'bse': {}, 'bands': {}, 'screening': {}}
    for chain in M.CHAINS:
        data['qp'][chain] = {lv: qp_payload(chain, lv) for lv in ('G0W0', 'scGW') if (chain, lv) not in WITHHELD['qp']}
        data['bse'][chain] = {run: bse_payload(chain, run) for run in ('BSE_G0W0', 'BSE_scGW')
                              if (chain, run) not in WITHHELD['bse']}
        data['bands'][chain] = {'soc': bands_payload(chain, True), 'nosoc': bands_payload(chain, False)}
        data['screening'][chain] = screening_payload(chain)
    # consistency assertions on the archive
    for chain in M.CHAINS:
        rows = M.read_qp(chain, 'G0W0')['rows']
        b = M.read_bands(chain, True)
        kg = np.arange(M.NK) / 48
        g_soc = np.interp(kg, b['k_frac'], b['energies'][M.VBM[chain]] - b['energies'][M.VBM[chain] - 1])
        g_eo = M.band(rows, M.VBM[chain] + 1, qp=False) - M.band(rows, M.VBM[chain], qp=False)
        assert np.abs(g_soc - g_eo).max() < 0.006, chain           # GW starting point = SOC band structure
        for run in ('BSE_G0W0', 'BSE_scGW'):
            assert len(M.read_bse(chain, run)['data']) == 2000
    for chain in ('Te3', 'As2S3', 'Bi2Te3'):                        # identical IP spectra: same QP energies
        a, s = M.read_bse(chain, 'BSE_G0W0')['data'], M.read_bse(chain, 'BSE_scGW')['data']
        assert np.abs(a[:, 3] - s[:, 3]).max() < 1e-5 * s[:, 3].max(), chain
        assert any('KfnQPdb' in l and M.SC_FILE[chain] in l for l in M.read_bse(chain, 'BSE_G0W0')['header'])
    assert any('KfnQPdb' in l and 'G0W0' in l for l in M.read_bse('S3', 'BSE_G0W0')['header'])

    # ------------------------------------------------------------ reference values
    gap = {c: {lv: M.qp_direct_gap(c, lv)[0] for lv in ('G0W0', 'scGW')} for c in M.CHAINS}
    dsoc = {c: M.soc_gap_change(c) for c in M.CHAINS}
    fb = {c: {run: M.first_bright(c, run) for run in ('BSE_G0W0', 'BSE_scGW')} for c in M.CHAINS}
    ref_gaps = {c: {'with_soc': gap[c]['scGW'], 'without_soc': gap[c]['scGW'] - dsoc[c]} for c in ('Te3', 'Bi2Te3')}
    ref_optical = {c: fb[c]['BSE_scGW'] for c in ('S3', 'Te3', 'As2S3')}
    sc_shift = {c: gap[c]['scGW'] - gap[c]['G0W0'] for c in M.CHAINS}
    ref_dw = {'S3': fb['S3']['BSE_scGW'] - fb['S3']['BSE_G0W0'] - sc_shift['S3'],
              'Te3': fb['Te3']['BSE_scGW'] - fb['Te3']['BSE_G0W0'],
              'As2S3': fb['As2S3']['BSE_scGW'] - fb['As2S3']['BSE_G0W0']}
    ref_eb = gap['Bi2Te3']['scGW'] - fb['Bi2Te3']['BSE_scGW']

    # ------------------------------------------------------------ wrong paths (ablations)
    mean_inc = float(np.mean([sc_shift[c] for c in ('S3', 'Te3', 'As2S3')]))
    eb_released = {c: gap[c]['scGW'] - fb[c]['BSE_scGW'] for c in ('S3', 'Te3', 'As2S3')}
    def first_feature(chain, run, frac=0.01):
        return M.first_bright(chain, run, fraction=frac)
    ff = {c: {run: first_feature(c, run) for run in ('BSE_G0W0', 'BSE_scGW')} for c in ('S3', 'As2S3')}
    k_bi = M.qp_direct_gap('Bi2Te3', 'G0W0')[1]
    me, mh = M.effective_masses('Bi2Te3', 'G0W0', k_bi)
    wire_bi = M.wire_exciton_binding(1 / (1 / me + 1 / mh), M.wire_polarizability('Bi2Te3'), M.quantum_radius('Bi2Te3'))
    calib = []
    wire_model = {}
    for c in ('S3', 'Te3'):
        m1, m2 = M.effective_masses(c, 'scGW', M.qp_direct_gap(c, 'scGW')[1])
        wire_model[c] = M.wire_exciton_binding(1 / (1 / m1 + 1 / m2), M.wire_polarizability(c), M.quantum_radius(c))
        calib.append(eb_released[c] / wire_model[c])
    abl = {
        'labels_trusted_dw': {'S3': ref_dw['S3'], 'Te3': ref_dw['Te3'] - sc_shift['Te3'],
                              'As2S3': ref_dw['As2S3'] - sc_shift['As2S3']},
        'labels_trusted_eb': gap['Bi2Te3']['G0W0'] - fb['Bi2Te3']['BSE_G0W0'],
        's3_same_energies_dw': fb['S3']['BSE_scGW'] - fb['S3']['BSE_G0W0'],
        'soc_added_again_gaps': {c: {'with_soc': gap[c]['scGW'] + dsoc[c], 'without_soc': gap[c]['scGW']}
                                 for c in ('Te3', 'Bi2Te3')},
        'increment_transfer_gap': gap['Bi2Te3']['G0W0'] + mean_inc,
        'increment_transfer_eb': gap['Bi2Te3']['G0W0'] + mean_inc - (fb['Bi2Te3']['BSE_G0W0'] + mean_inc),
        'g0w0_gap_for_bi2te3': gap['Bi2Te3']['G0W0'],
        'first_feature_optical_As2S3': ff['As2S3']['BSE_scGW'],
        'first_feature_dw': {'S3': ff['S3']['BSE_scGW'] - ff['S3']['BSE_G0W0'] - sc_shift['S3'],
                             'Te3': ref_dw['Te3'], 'As2S3': ff['As2S3']['BSE_scGW'] - ff['As2S3']['BSE_G0W0']},
        'wire_model_eb': wire_bi,
        'wire_model_calibrated_eb': wire_bi * float(np.mean(calib)),
        'released_mean_eb': float(np.mean(list(eb_released.values()))),
    }
    # defensible variations (robustness)
    ip_gap_peak = M.ip_first_peak('Bi2Te3', 'BSE_G0W0')
    ip_gap_onset = M.ip_onset('Bi2Te3', 'BSE_G0W0')
    s3_ip_shift = M.ip_first_peak('S3', 'BSE_scGW') - M.ip_first_peak('S3', 'BSE_G0W0')
    eb_w0 = gap['Bi2Te3']['scGW'] - fb['Bi2Te3']['BSE_G0W0']
    eb_w0_ip = ip_gap_peak - fb['Bi2Te3']['BSE_G0W0']
    ratios = {c: (gap[c]['scGW'] - ref_optical[c]) / (gap[c]['scGW'] - ref_optical[c] + ref_dw[c]) for c in ref_dw}
    rel_sc = {c: sc_shift[c] / gap[c]['G0W0'] for c in ref_dw}
    slope = float(np.mean([(ratios[c] - 1) / rel_sc[c] for c in ref_dw]))
    rel_bi = (ip_gap_peak - gap['Bi2Te3']['G0W0']) / gap['Bi2Te3']['G0W0']
    # Te3 SOC-free gap from k-resolved corrections applied to the no-SOC bands
    rows = M.read_qp('Te3', 'scGW')['rows']
    dq = (M.band(rows, 49) - M.band(rows, 48)) - (M.band(rows, 49, False) - M.band(rows, 48, False))
    nb = M.read_bands('Te3', False)
    g_nosoc = np.interp(np.arange(M.NK) / 48, nb['k_frac'], nb['energies'][M.VBM_NOSOC['Te3']] - nb['energies'][M.VBM_NOSOC['Te3'] - 1])
    def soc_change_at_k(chain, kidx):
        kf = (kidx - 1) / 48
        out = []
        for soc, v in ((True, M.VBM[chain]), (False, M.VBM_NOSOC[chain])):
            b = M.read_bands(chain, soc)
            out.append(float(np.interp(kf, b['k_frac'], b['energies'][v] - b['energies'][v - 1])))
        return out[0] - out[1]
    te3_rows = M.read_qp('Te3', 'scGW')['rows']
    te3_fund = float(M.band(te3_rows, 49).min() - M.band(te3_rows, 48).max())
    rob = {
        'te3_without_soc_fixed_k': gap['Te3']['scGW'] - soc_change_at_k('Te3', M.qp_direct_gap('Te3', 'scGW')[1]),
        'bi2te3_without_soc_fixed_k': gap['Bi2Te3']['scGW'] - soc_change_at_k('Bi2Te3', M.qp_direct_gap('Bi2Te3', 'scGW')[1]),
        'te3_fundamental_gap': te3_fund,
        'bi2te3_gap_ip_first_peak': ip_gap_peak, 'bi2te3_gap_ip_onset': ip_gap_onset,
        's3_dw_ip_shift': fb['S3']['BSE_scGW'] - fb['S3']['BSE_G0W0'] - s3_ip_shift,
        'te3_without_soc_k_resolved': float((g_nosoc + dq).min()),
        'eb_no_w_correction': eb_w0_ip,
        'eb_additive_mean_dw': ip_gap_peak - (fb['Bi2Te3']['BSE_G0W0'] + float(np.mean(list(ref_dw.values())))),
        'eb_mean_ratio': eb_w0_ip * float(np.mean(list(ratios.values()))),
        'eb_ratio_scaled_by_self_consistency': eb_w0_ip * (1 + slope * rel_bi),
    }
    hidden = {'targets': {'scgw_direct_gap_eV': ref_gaps, 'scgw_optical_gap_eV': ref_optical,
                          'screening_update_shift_eV': ref_dw, 'bi2te3_binding_energy_eV': ref_eb},
              'thresholds': THRESHOLDS,
              'reference_notes': {'bi2te3_scgw_gap_k_index': M.qp_direct_gap('Bi2Te3', 'scGW')[1],
                                  'bi2te3_scgw_optical_gap_eV': fb['Bi2Te3']['BSE_scGW'],
                                  'soc_gap_change_eV': dsoc, 'self_consistency_gap_shift_eV': sc_shift}}
    values = [ref_gaps['Te3']['with_soc'], ref_gaps['Te3']['without_soc'], ref_gaps['Bi2Te3']['with_soc'],
              ref_gaps['Bi2Te3']['without_soc'], *ref_optical.values(), *ref_dw.values(), ref_eb,
              fb['Bi2Te3']['BSE_scGW'], sc_shift['Bi2Te3']]

    catalog = {
        'title': 'Exfoliable 1D chalcogenide chains: archived DFT, GW and BSE results (Grillo et al., ACS Nano 2026)',
        'chains': {c: {'atoms_per_cell': FORMULA_ATOMS[c], 'chain_lattice_constant_A': M.CHAIN_C_A[c],
                       'valence_band_max_index': M.VBM[c],
                       'qp_band_range_archived': [int(M.read_qp(c, 'G0W0')['rows'][:, 1].min()),
                                                  int(M.read_qp(c, 'G0W0')['rows'][:, 1].max())],
                       'scGW_iteration': M.SC_FILE[c]} for c in M.CHAINS},
        'k_points': {'count': M.NK, 'k_frac_of_index': {str(n): round((n - 1) / 48, 6) for n in range(1, M.NK + 1)},
                     'note': 'irreducible points of a 48-point grid along the chain; k_frac in units of 2*pi/c'},
        'levels': {'G0W0': 'one-shot GW on PBE (archived o-G0W0.qp)',
                   'scGW': 'last archived iteration of eigenvalue self-consistent GW, G and W updated (o-G3W3.qp; o-G4W4.qp for As2S3)'},
        'requests': {
            'qp': {'example': {'kind': 'qp', 'chain': 'Te3', 'level': 'scGW', 'k_indices': [20, 21, 22]},
                   'cost': 'one unit per k-point', 'returns': 'all archived bands at those k-points (Eo, E-Eo, Sc|Eo) and the Yambo run header'},
            'bse_spectrum': {'example': {'kind': 'bse_spectrum', 'chain': 'Te3', 'run': 'BSE_scGW'},
                             'runs': ['BSE_G0W0', 'BSE_scGW'], 'cost': COSTS['bse_spectrum'],
                             'returns': 'the archived Haydock spectrum (0-10 eV, field along the chain) with its Yambo header; run names are the authors\' labels'},
            'dft_bands': {'example': {'kind': 'dft_bands', 'chain': 'Te3', 'soc': True}, 'cost': COSTS['dft_bands'],
                          'returns': 'archived PBE bands along Gamma-Z (301 points), computed with or without spin-orbit coupling, bands near the gap'},
            'wire_screening': {'example': {'kind': 'wire_screening', 'chain': 'Te3'}, 'cost': COSTS['wire_screening'],
                               'returns': 'DFPT electronic dielectric tensor of the supercell, cells and the Environ quantum-volume record'}},
        'withheld_for_grading': ['qp of Bi2Te3 at the scGW level', 'bse_spectrum Bi2Te3 BSE_scGW'],
        'purchase_costs': COSTS, 'budget': BUDGET,
        'methods': ('PBE with Quantum ESPRESSO; GW and BSE with Yambo 5.2.3 (plasmon-pole GW, cylindrical Coulomb cutoff, '
                    'random-integration method; BSE in the Tamm-Dancoff approximation solved by the Haydock method with '
                    '5 meV damping, light polarized along the chain). Isolated chains in periodic supercells with vacuum '
                    'perpendicular to the chain. The headers returned with each purchase record the run parameters.'),
    }
    (TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')
    (TASK/'environment/lab/data.json').write_text(json.dumps(data, separators=(',', ':')) + '\n')
    (TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
    (TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(values) + '\n')
    near_miss = {'te3_gap_from_first_ip_peak': M.ip_first_peak('Te3', 'BSE_scGW')}
    (TASK/'solution/ablation_values.json').write_text(json.dumps({'ablations': abl, 'robustness': rob, 'near_miss': near_miss},
                                                                 indent=1) + '\n')
    for path, extra in ((TASK/'environment/lab/public.json', {'reserved_validation_charge': 0}),
                        (TASK/'tests/public.json', {'min_report_chars': 400})):
        doc = json.loads(path.read_text())
        doc.update(solver_budget=BUDGET, **extra)
        path.write_text(json.dumps(doc, indent=1) + '\n')
    summary = {'reference': hidden['targets'], 'notes': hidden['reference_notes'], 'ablations': abl, 'robustness': rob,
               'wire_model_released': wire_model, 'calibration_factors': calib, 'ratios': ratios, 'slope': slope,
               'catalog_units': sum(COSTS['qp_per_k_point'] * M.NK for c in M.CHAINS for lv in ('G0W0', 'scGW')
                                    if (c, lv) not in WITHHELD['qp']) + 7 * COSTS['bse_spectrum'] + 8 + 4,
               'data_json_MB': (TASK/'environment/lab/data.json').stat().st_size / 2**20}
    print(json.dumps(summary, indent=1, default=float))


if __name__ == '__main__':
    main()
