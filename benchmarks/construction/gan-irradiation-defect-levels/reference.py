"""Reference analysis and wrong-path variants for gan-irradiation-defect-levels (numpy only).

`analyse(**options)` returns the graded quantities for one analysis protocol. The defaults are the reference
protocol; each option reproduces one defensible variation or one plausible wrong path.

Reference protocol
- Host energy: the perfect supercell computed like the defect cells (k = (1/2,1/2,1/2), Ga PAW set).
- Band edges: Gamma-point eigenvalues of the Gamma-only perfect supercell with the same PAW set
  (perfect/kpt_gam); its eigenvalue scale is used as is.
- Finite size: eFNV (Kumagai-Oba) with the static (relaxed-ion) DFPT tensor, PC + site-potential alignment
  over atoms outside the sphere inscribed in the supercell around the defect.
- Configurations: for each interstitial charge state the lower of the two archived sites.
- Levels: thermodynamic transitions between ground states for E_F between the VBM and the CBM.
"""
import json
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
POOL = M.load_pool()
DIELECTRIC = json.loads((HERE/'dielectric.json').read_text())
DEFECTS = {'V_Ga': ['Va_Ga1'], 'V_N': ['Va_N1'], 'Ga_i': ['Ga_i1', 'Ga_i2'], 'N_i': ['N_i1', 'N_i2'],
           'C_N': ['C_N1'], 'Ga_N': ['Ga_N1'], 'N_Ga': ['N_Ga1'], 'VV': ['VV_GaN1']}
GRADED_LEVEL_DEFECTS = ['V_Ga', 'V_N', 'Ga_i', 'N_i', 'C_N']
EF_OFFSET_BELOW_CBM = 0.5


def runs_of(base):
    return {M.split_name(n)[1]: n for n in POOL if not n.startswith('perfect') and M.split_name(n)[0] == base}


def band_edges(source='gamma'):
    """(VBM, CBM) on the perfect-supercell eigenvalue scale for several (right and wrong) choices."""
    def edges(name):
        ev = POOL[name]['eigenvalues']
        occ = [e for s in ev for e, o in s if o > 0.5]
        emp = [e for s in ev for e, o in s if o < 0.5]
        return max(occ), min(emp)
    if source == 'gamma':
        return edges('perfect/kpt_gam')
    if source == 'gamma_aligned_all_sites':
        v, c = edges('perfect/kpt_gam')
        shift = float(np.mean(np.array(POOL['perfect']['core_potentials_V']) - np.array(POOL['perfect/kpt_gam']['core_potentials_V'])))
        return v + shift, c + shift
    if source == 'gamma_aligned_N_sites':
        v, c = edges('perfect/kpt_gam')
        sym = np.array(POOL['perfect']['symbols'])
        d = np.array(POOL['perfect']['core_potentials_V']) - np.array(POOL['perfect/kpt_gam']['core_potentials_V'])
        return v + float(d[sym == 'N'].mean()), c + float(d[sym == 'N'].mean())
    if source == 'R_point':
        return edges('perfect')
    if source == 'efermi':
        return POOL['perfect']['efermi_eV'], edges('perfect/kpt_gam')[1]
    if source == 'gamma_Ga_d':
        return edges('perfect/WAV_ref')
    raise ValueError(source)


def corrected_energies(base, scheme='efnv', eps='eps_0', region='outside_inscribed', host='perfect',
                       exclude=(), ewald_cache={}):
    """{q: E_defect(q) + E_corr(q) - E_host} for one archived configuration family."""
    perfect = POOL['perfect']
    host_e = POOL[host]['energy_sigma0_eV']
    if eps == 'eps_0_iso':
        tensor = np.eye(3) * np.trace(np.array(DIELECTRIC['eps_0'])) / 3
    else:
        tensor = np.array(DIELECTRIC[eps]) if eps in DIELECTRIC else np.eye(3) * float(eps)
    key = (eps, tensor.tobytes())
    if key not in ewald_cache:
        ewald_cache[key] = M.Ewald(perfect['cell_A'], tensor)
    ew = ewald_cache[key]
    out = {}
    for q, name in runs_of(base).items():
        if name in exclude:
            continue
        rec = POOL[name]
        if scheme == 'none' or q == 0:
            corr = 0.0
        elif scheme == 'efnv':
            corr = M.efnv(rec, perfect, q, tensor, region=region, ewald=ew)['total']
        elif scheme == 'efnv_wrong_sign':
            corr = M.efnv(rec, perfect, q, tensor, region=region, ewald=ew, potential_sign=+1)['total']
        elif scheme == 'pc_only':
            corr = M.efnv(rec, perfect, q, tensor, alignment=False, ewald=ew)['total']
        else:
            raise ValueError(scheme)
        out[q] = rec['energy_sigma0_eV'] + corr - host_e
    return out


def defect_energies(defect, sites='min', **kw):
    fams = DEFECTS[defect]
    if sites != 'min':
        fams = [sites] if sites in fams else fams[:1]
    per = [corrected_energies(f, **kw) for f in fams]
    qs = sorted(set().union(*per))
    return {q: min(p[q] for p in per if q in p) for q in qs}


def analyse(edges='gamma', scheme='efnv', eps='eps_0', region='outside_inscribed', host='perfect',
            sites=None, exclude=(), binding_charge_states='ground', fp_charge_states='ground'):
    sites = sites or {}
    vbm, cbm = band_edges(edges)
    gap = cbm - vbm
    kw = dict(scheme=scheme, eps=eps, region=region, host=host, exclude=exclude)
    E = {d: defect_energies(d, sites=sites.get(d, 'min'), **kw) for d in DEFECTS}
    levels = {}
    for d in GRADED_LEVEL_DEFECTS + ['VV', 'Ga_N', 'N_Ga']:
        inside, _ = M.levels_from_energies(E[d], vbm, gap)
        levels[d] = [{'from': a, 'to': b, 'level_eV': x} for (a, b), x in inside]
    ef = gap - EF_OFFSET_BELOW_CBM

    def ef_state(d):
        if binding_charge_states == 'neutral':
            return 0, E[d][0]
        return M.ground_state(E[d], vbm, ef)
    (q1, f1), (q2, f2), (q3, f3) = ef_state('V_Ga'), ef_state('V_N'), ef_state('VV')
    binding = f1 + f2 - f3

    def fp(dv, di):
        if fp_charge_states == 'neutral':
            return E[dv][0] + E[di][0], (0, 0)
        (a, fa), (b, fb) = M.ground_state(E[dv], vbm, ef), M.ground_state(E[di], vbm, ef)
        return fa + fb, (a, b)
    fpn, qn = fp('V_N', 'N_i')
    fpg, qg = fp('V_Ga', 'Ga_i')
    return {'vbm': vbm, 'cbm': cbm, 'gap_eV': gap, 'levels': levels,
            'binding_energy_eV': binding, 'binding_states': [q1, q2, q3],
            'frenkel_pair_energy_eV': {'N': fpn, 'Ga': fpg}, 'frenkel_states': {'N': qn, 'Ga': qg},
            'energies': {d: {str(q): v for q, v in e.items()} for d, e in E.items()}}


if __name__ == '__main__':
    r = analyse()
    print(json.dumps({k: v for k, v in r.items() if k != 'energies'}, indent=1))
