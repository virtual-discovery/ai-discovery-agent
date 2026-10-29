"""Build the staged, lab, private and ablation files of hgte-gated-subband-densities.

    bash run_in_image.sh make_task.py      (needs /scratch/pool.npz from extract.py and the archive at /data)

Writes into benchmarks/tasks/hgte-gated-subband-densities/:
  environment/agent/catalog.json         free agent-visible metadata
  environment/lab/data/lab_pool.npz      transport grids, Hall scans, trimmed kdotpy runs (served by the lab)
  environment/lab/data/calcs.json        run metadata (archived .kpout text, folders)
  environment/lab/private/hidden.json    graded targets and tolerances
  environment/lab/private/hidden_values.json
  solution/ablation_values.json          wrong paths and defensible variants
and calibration.json here.
"""
import json
from pathlib import Path

import numpy as np

import sdh
from model import Pool, interp

HERE = Path(__file__).resolve().parent
TASK = HERE.parents[1]/'tasks/hgte-gated-subband-densities'
ARCH = Path('/data/calc/data_repository')
P = Pool()

SAMPLES = {'26nm': '26', '70nm': '70', '107nm': '107'}
RPRE = {'26': 1000.0, '70': 557.0, '107': 10000.0}
EWIN = (-60.0, 80.0)          # energy window (meV) of the integrated DOS served by the lab
GATE_WINDOW = (0.05, 0.29)    # electron-side window of the gate calibration (authors' 70 nm convention)
HALL_BMAX = 0.2

# ---------------------------------------------------------------- gate calibration (70 nm)
tr = {s: (P.P[f'tr_{s}_vg'], P.P[f'tr_{s}_B'], P.P[f'tr_{s}_data']) for s in SAMPLES.values()}


def hall_line(s, bmax=HALL_BMAX, window=GATE_WINDOW):
    vg, B, D = tr[s]
    nH = np.array([sdh.hall_density(B, D[g, :, 0], D[g, :, 2], RPRE[s], bmax) for g in range(len(vg))])
    m = (vg > window[0] + 1e-9) & (vg <= window[1] + 1e-9)
    a, b = np.polyfit(vg[m], nH[m], 1)
    return float(a), float(b), nH


A70, B70, NH70 = hall_line('70')
assert abs(A70 - 9.81968542) < 1e-6 and abs(B70 - 3.20604992) < 1e-6     # authors' FigureS3 values
A70_b01, B70_b01, _ = hall_line('70', bmax=0.1)


def n_of_vg(v, a=A70, b=B70):
    return a * v + b          # 1e11 cm^-2


# ---------------------------------------------------------------- model (archived kdotpy runs, 70 nm)
NM = P.n_values('70', 'main')
NW = P.n_values('70', 'w165')
BD = {n: P.band_densities('70', n * 1) for n in NM}


def model_band(band, v, a=A70, b=B70):
    return interp(n_of_vg(v, a, b), [n * 1e3 for n in NM], [BD[n][band] for n in NM])


def model_pocket(v, wedge='w165', own=False, a=A70, b=B70):
    return interp(n_of_vg(v, a, b), [n * 1e3 for n in NW], [P.pocket('70', n, wedge, own) for n in NW])


M1 = model_band(2, 0.20)          # bottom surface (band 2) on the electron side
M3 = model_band(2, -0.45)         # top surface: band 2 on the hole side (labels cross near n ~ 0.001 nm^-2)
M2 = model_pocket(-0.60)          # majority ([11]-type) pocket, single pocket
# consistency: main-run band -1 = 2 x (majority + minority pocket)
for n in NW:
    full = -BD[n][-1]
    two = 2 * (P.pocket('70', n, 'w165') + P.pocket('70', n, 'w105'))
    assert full < 0.5 or abs(two / full - 1) < 0.02, (n, full, two)

# ---------------------------------------------------------------- experiment (authors' 70 nm chart)
a = np.genfromtxt(ARCH/'FFT_experiment/70nm_FT.csv', delimiter=',', skip_header=1)[:, 1:]
VG = np.unique(a[:, 0])
NN = np.unique(a[:, 1])
PW = np.full((len(VG), len(NN)), np.nan)
PW[np.searchsorted(VG, a[:, 0]), np.searchsorted(NN, a[:, 1])] = a[:, 2]
assert not np.isnan(PW).any()
HW = 0.02
# (gate voltage, ridge guess from the hand-tracked 2D chart, search window); see README "Ridge assignment"
RIDGES = {'top_surface_Vg_p0.20': (0.20, 2.81, 0.3), 'bottom_surface_Vg_p0.05': (0.05, 0.83, 0.15),
          'bulk_subband_E1_Vg_p0.30': (0.30, 1.22, 0.15), 'top_surface_Vg_m0.45': (-0.45, 0.45, 0.1),
          'bottom_surface_Vg_m0.45': (-0.45, 0.63, 0.1), 'hole_pocket_Vg_m0.70': (-0.70, 2.49, 0.2)}
HOLE_SLOPE, HOLE_V0 = 5.2, -0.22    # hole ridge n = 5.2 (V0 - Vg), from a line search on the chart


def ridge_values(vg, n, pw, hw):
    out = {}
    for k, (v, g, w) in RIDGES.items():
        out[k] = sdh.peak_near(n, sdh.averaged_log(vg, pw, v, hw), g, w)
    pts = []
    for v in np.arange(-0.75, -0.549, 0.05):
        pk = sdh.peak_near(n, sdh.averaged_log(vg, pw, v, hw), HOLE_SLOPE * (HOLE_V0 - v), 0.2)
        pts.append((v, pk))
    sl, ic = np.polyfit([p[0] for p in pts], [p[1] for p in pts], 1)
    out['hole_onset_Vg'] = -ic / sl
    return out, pts


EXP, HOLE_PTS = ridge_values(VG, NN, PW, HW)


def variant(quantity='sxx', hw=HW, **kw):
    vg, B, D = tr['70']
    n, pw = sdh.chart(vg, B, D, RPRE['70'], quantity=quantity, **kw)
    return ridge_values(vg, n, pw, hw)[0]


# exact reproduction check of the archived chart from the raw record
n_chk, pw_chk = sdh.chart(*tr['70'][:2], tr['70'][2], RPRE['70'])
assert np.allclose(n_chk, NN) and np.max(np.abs(pw_chk - PW)) / PW.max() < 1e-10

VARIANTS = {'single_spectrum': variant(hw=0.0), 'bmax1p0_padded': variant(bmax=1.0, pad=4),
            'rho_xx_hann': variant(quantity='rxx', window='hann'), 'block_average': variant(hw=0.024),
            'bmax0p6_padded': variant(bmax=0.6, pad=4), 'bmax1p5_bmin0p2_padded': variant(bmax=1.5, bmin=0.2, pad=4)}

# ---------------------------------------------------------------- ablations (wrong paths from the data)
L70 = sdh.averaged_log(VG, PW, -0.70, HW)
maxima = [p for p in sdh.local_maxima(NN, L70) if 1.0 <= p[0] <= 3.2]
hole_strongest = max(maxima, key=lambda p: p[1])[0]
model_hole_m070 = model_pocket(-0.70)
hole_nearest_model = min(maxima, key=lambda p: abs(p[0] - model_hole_m070))[0]
sign_change = float(VG[np.where(np.diff(np.sign(NH70)))[0][0]])
rigid_bd, rigid_E = P.rigid('70', 0.0, n_of_vg(0.20))
# positive-gate misassignment seen in the dev pilots: the steep top-surface ridge taken as the first bulk
# subband, and a flat feature near 2x the bottom surface (harmonic) or a sum frequency taken as the top surface
L20, L30 = sdh.averaged_log(VG, PW, 0.20, HW), sdh.averaged_log(VG, PW, 0.30, HW)
harm_p020 = sdh.peak_near(NN, L20, 1.5, 0.1)
sum_p020 = sdh.peak_near(NN, L20, 1.75, 0.1)
steep_p030 = sdh.peak_near(NN, L30, 3.17, 0.15)
ABL = {
    'steep_ridge_as_bulk_harmonic_as_top': {'top_surface_Vg_p0.20': harm_p020, 'bulk_subband_E1_Vg_p0.30': steep_p030},
    'steep_ridge_as_bulk_sum_as_top': {'top_surface_Vg_p0.20': sum_p020, 'bulk_subband_E1_Vg_p0.30': steep_p030},
    'bulk_spin_degenerate': {'bulk_subband_E1_Vg_p0.30': 2 * EXP['bulk_subband_E1_Vg_p0.30'],
                             'hole_pocket_Vg_m0.70': 2 * EXP['hole_pocket_Vg_m0.70']},
    'hole_strongest_peak': {'hole_pocket_Vg_m0.70': hole_strongest},
    'hole_peak_nearest_model': {'hole_pocket_Vg_m0.70': hole_nearest_model},
    'top_is_larger_density': {'top_surface_Vg_m0.45': EXP['bottom_surface_Vg_m0.45'],
                              'bottom_surface_Vg_m0.45': EXP['top_surface_Vg_m0.45'],
                              'model_top_surface_Vg_m0.45': model_band(1, -0.45)},
    'four_equivalent_pockets': {'majority_pocket_multiplicity': 4},
    'single_hole_pocket': {'majority_pocket_multiplicity': 1},
    'wedge_own_fermi_level': {'model_hole_pocket_Vg_m0.60': model_pocket(-0.60, own=True)},
    'band_minus1_split_in_two': {'model_hole_pocket_Vg_m0.60': -model_band(-1, -0.60) / 2},
    'band_minus1_split_in_four': {'model_hole_pocket_Vg_m0.60': -model_band(-1, -0.60) / 4},
    'wedge_not_normalised': {'model_hole_pocket_Vg_m0.60': model_pocket(-0.60) * 3},
    'rigid_band_from_neutral_run': {'model_bottom_surface_Vg_p0.20': rigid_bd[2]},
    'holes_appear_at_neutrality': {'hole_onset_Vg': -B70 / A70},
    'holes_appear_at_hall_sign_change': {'hole_onset_Vg': sign_change},
}

# ---------------------------------------------------------------- targets and tolerances
TARGETS = dict(EXP)
TARGETS.update({'model_bottom_surface_Vg_p0.20': M1, 'model_top_surface_Vg_m0.45': M3,
                'model_hole_pocket_Vg_m0.60': M2})
TOL = {k: {'rel': 0.06, 'abs': 0.07} for k in RIDGES}
TOL['hole_onset_Vg'] = {'rel': 0.0, 'abs': 0.06}
for k in ('model_bottom_surface_Vg_p0.20', 'model_top_surface_Vg_m0.45', 'model_hole_pocket_Vg_m0.60'):
    TOL[k] = {'rel': 0.06, 'abs': 0.04}


def passes(key, x):
    t = TARGETS[key]
    return abs(x - t) <= max(TOL[key]['rel'] * abs(t), TOL[key]['abs'])


ROB = {name: dict(v) for name, v in VARIANTS.items()}
A70_w, B70_w, _ = hall_line('70', window=(0.06, 0.28))
ROB['hall_fit_window_0p06_0p28'] = {'model_bottom_surface_Vg_p0.20': model_band(2, 0.20, A70_w, B70_w),
                                    'model_top_surface_Vg_m0.45': model_band(2, -0.45, A70_w, B70_w),
                                    'model_hole_pocket_Vg_m0.60': model_pocket(-0.60, a=A70_w, b=B70_w)}


def quad(x, xs, ys):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    i = np.argsort(np.abs(xs - x))[:3]
    return float(np.polyval(np.polyfit(xs[i], ys[i], 2), x))


ROB['quadratic_interpolation'] = {
    'model_bottom_surface_Vg_p0.20': quad(n_of_vg(0.20), [n * 1e3 for n in NM], [BD[n][2] for n in NM]),
    'model_top_surface_Vg_m0.45': quad(n_of_vg(-0.45), [n * 1e3 for n in NM], [BD[n][2] for n in NM]),
    'model_hole_pocket_Vg_m0.60': quad(n_of_vg(-0.60), [n * 1e3 for n in NW], [P.pocket('70', n, 'w165') for n in NW])}
# model hole onset: majority pocket density extrapolated to zero
pk = [(n * 1e3, P.pocket('70', n, 'w165')) for n in NW if P.pocket('70', n, 'w165') > 0.01]
sl, ic = np.polyfit([p[0] for p in pk[-3:]], [p[1] for p in pk[-3:]], 1)
model_onset = (-ic / sl - B70) / A70
ROB['model_pocket_onset'] = {'hole_onset_Vg': model_onset}

PROBLEMS = []
for name, vals in ROB.items():
    PROBLEMS += [('robustness fails', name, k, x, TARGETS[k]) for k, x in vals.items() if not passes(k, x)]
for name, vals in ABL.items():
    PROBLEMS += [('ablation passes', name, k, x, TARGETS[k]) for k, x in vals.items()
                 if k != 'majority_pocket_multiplicity' and passes(k, x)]

# ---------------------------------------------------------------- lab data
data_dir = TASK/'environment/lab/data'
data_dir.mkdir(parents=True, exist_ok=True)
arrays, calcs = {}, {}
for name, s in SAMPLES.items():
    vg, B, D = tr[s]
    arrays[f'tr_{name}_vg'], arrays[f'tr_{name}_B'], arrays[f'tr_{name}_data'] = vg, B, D
    m = np.abs(B) <= HALL_BMAX + 1e-9
    slopes = np.array([np.polyfit(B[m], D[g, m, 2] / D[g, m, 0], 1)[0] for g in range(len(vg))])
    i0 = int(np.abs(B).argmin())
    arrays[f'hall_{name}'] = np.c_[vg, slopes, D[:, i0, 1] / D[:, i0, 0]]
for key, info in P.meta['calcs'].items():
    E, I = P.P[f'E_{key}'], P.P[f'I_{key}']
    w = (E >= EWIN[0] - 1e-9) & (E <= EWIN[1] + 1e-9)
    name = [k for k, v in SAMPLES.items() if v == info['sample']][0]
    kind = {'main': 'main', 'w105': 'kphi_105_165', 'w165': 'kphi_165_285'}[info['kind']]
    lk = f'{name}|{kind}|{info["n_nm2"]:+.4f}'
    arrays[f'E|{lk}'], arrays[f'I|{lk}'], arrays[f'P|{lk}'] = E[w], I[w], P.P[f'P_{key}']
    calcs[lk] = {'sample': name, 'run': kind, 'n_nm2': info['n_nm2'], 'kpout': info['kpout_text'],
                 'archive_folder': 'FigureS3/calculations/' + info['folder']}
np.savez_compressed(data_dir/'lab_pool.npz', **arrays)
(data_dir/'calcs.json').write_text(json.dumps(calcs, indent=1) + '\n')

# ---------------------------------------------------------------- agent catalog
centers = {}
for name, s in SAMPLES.items():
    vg = tr[s][0]
    centers[name] = [round(c, 2) for c in np.arange(-0.8, 0.651, 0.05) if vg.min() - 1e-9 <= c <= vg.max() + 1e-9]
readme = (ARCH/'FigureS3/README.md').read_text()
catalog = {
    'description': 'Gated HgTe/Hg0.32Cd0.68Te topological-insulator films on CdTe(001): raw mK magnetotransport of three '
                   'Hall bars and the authors\' released self-consistent 8-band k.p (kdotpy) runs for the same layer '
                   'thicknesses. Everything listed here can be purchased from the lab; nothing else is available.',
    'samples': {
        '26nm': {'HgTe_thickness_nm': 26, 'transport_record': '26nm.csv', 'pre_resistor_ohm': 1000,
                 'top_barrier_HgCdTe_nm': 50, 'HfOx_gate_dielectric_nm': 13.6, 'cryostat': '"18T", base ~17 mK'},
        '70nm': {'HgTe_thickness_nm': 70, 'transport_record': '70nm.csv', 'pre_resistor_ohm': 557,
                 'top_barrier_HgCdTe_nm': 45, 'HfOx_gate_dielectric_nm': 15, 'cryostat': '"18T", base ~17 mK'},
        '107nm': {'HgTe_thickness_nm': 107, 'transport_record': '110nm.csv (nominal 110 nm; modelled as 107 nm)',
                  'pre_resistor_ohm': 10000, 'top_barrier_HgCdTe_nm': 44, 'HfOx_gate_dielectric_nm': 14.2,
                  'cryostat': '"18T", base ~17 mK'},
    },
    'device': 'Top-gated (HfOx / metal) Hall bars, 600 um x 200 um; the longitudinal contacts are 600 um apart '
              '(L/W = 3). Lock-in at 13.7 Hz, 200 uV excitation, in series with the listed pre-resistor; VI is the '
              'voltage across the pre-resistor. Field and voltage offsets are corrected in the records.',
    'transport': {'columns': ['VG (V)', 'B (T)', 'VI (V)', 'Vxx (V)', 'Vxy (V)'], 'field_range_T': 'record minimum to 2.0',
                  'gate_step_V': 0.002, 'block_half_width_V': 0.025, 'vg_centers': centers},
    'hall_scan': 'per gate voltage: slope of Vxy/VI versus B for |B| <= 0.2 T (1/T) and Vxx/VI at B ~ 0',
    'calculations': {name: P.n_values(s, 'main') for name, s in SAMPLES.items()},
    'pocket_calculations': {'n_nm2': {name: P.n_values(s, 'w165') for name, s in SAMPLES.items()},
                            'kphi_deg': [[105, 165], [165, 285]]},
    'calculation_output': 'archived .kpout text, potential profile V(z) (meV), and the band-resolved integrated DOS '
                          'IDOS(b) and IDOS(total) (1/nm^2) on the archived 0.1 meV grid for %g..%g meV' % EWIN,
    'kdotpy_readme_FigureS3': readme,
}
(TASK/'environment/agent/catalog.json').write_text(json.dumps(catalog, indent=1) + '\n')

# ---------------------------------------------------------------- private and ablation files
FIELDS = {'top_surface_Vg_p0.20': ['experiment_1e11_cm2', 'top_surface_at_+0.20V'],
          'bottom_surface_Vg_p0.05': ['experiment_1e11_cm2', 'bottom_surface_at_+0.05V'],
          'bulk_subband_E1_Vg_p0.30': ['experiment_1e11_cm2', 'first_bulk_subband_at_+0.30V'],
          'top_surface_Vg_m0.45': ['experiment_1e11_cm2', 'top_surface_at_-0.45V'],
          'bottom_surface_Vg_m0.45': ['experiment_1e11_cm2', 'bottom_surface_at_-0.45V'],
          'hole_pocket_Vg_m0.70': ['experiment_1e11_cm2', 'hole_pocket_at_-0.70V'],
          'hole_onset_Vg': [None, 'hole_onset_gate_V'],
          'model_bottom_surface_Vg_p0.20': ['model_1e11_cm2', 'bottom_surface_at_+0.20V'],
          'model_top_surface_Vg_m0.45': ['model_1e11_cm2', 'top_surface_at_-0.45V'],
          'model_hole_pocket_Vg_m0.60': ['model_1e11_cm2', 'hole_pocket_at_-0.60V'],
          'majority_pocket_multiplicity': [None, 'majority_pocket_multiplicity']}
hidden = {'fields': FIELDS, 'targets': TARGETS, 'tolerances': TOL, 'majority_pocket_multiplicity': 2,
          'reference': {'gate_line_70nm': [A70, B70], 'hole_ridge_points': HOLE_PTS}}
(TASK/'environment/lab/private/hidden.json').write_text(json.dumps(hidden, indent=1) + '\n')
(TASK/'environment/lab/private/hidden_values.json').write_text(json.dumps(list(TARGETS.values())) + '\n')
(TASK/'solution/ablation_values.json').write_text(json.dumps({'fields': FIELDS, 'ablations': ABL, 'robustness': ROB},
                                                              indent=1) + '\n')

calib = {}
for k, t in TARGETS.items():
    tol = max(TOL[k]['rel'] * abs(t), TOL[k]['abs'])
    calib[k] = {'reference': round(t, 4), 'tolerance': round(tol, 4),
                'robustness': {n: round(v[k], 4) for n, v in ROB.items() if k in v},
                'ablations': {n: round(v[k], 4) for n, v in ABL.items() if k in v}}
calib['majority_pocket_multiplicity'] = {'reference': 2, 'tolerance': 'exact',
                                         'ablations': {'four_equivalent_pockets': 4, 'single_hole_pocket': 1}}
calib['_notes'] = {'gate_line_70nm': [A70, B70], 'gate_line_B0p1_not_the_convention': [A70_b01, B70_b01], 'gate_line_window_0p06_0p28': [A70_w, B70_w], 'rigid_band_energy_meV': float(rigid_E),
                   'model_hole_pocket_Vg_m0.70': model_hole_m070, 'hole_ridge_points': HOLE_PTS,
                   'model_top_band1_m045': model_band(1, -0.45)}
(HERE/'calibration.json').write_text(json.dumps(calib, indent=1) + '\n')
print(json.dumps(calib, indent=1))
assert not PROBLEMS, PROBLEMS
