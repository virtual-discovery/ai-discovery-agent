"""Extract the archived EPW transport results and force constants into pool.json.

    bash run_in_image.sh extract.py /data        (DEST of fetch.py mounted at /data)

pool.json holds, unmodified apart from parsing:
  runs  - every archived EPW ibte calculation: input parameters (ibte.in), lattice parameter (scf.in), the fine
          meshes actually used (from the output), and for each temperature the Fermi level, the printed DOS, the
          SERTA conductivity tensor, the xx element and max error of every IBTE iteration, and the final IBTE tensor.
  ifcs  - every q2r (DFPT) and TDEP (converted to QE form) force-constant set: header geometry, mass, mesh and the
          force constants (Ry/bohr^2), plus for TDEP sets the sampling temperature and TDEP unit cell.
"""
import json
import re
import sys
from pathlib import Path

SRC = Path(sys.argv[1] if len(sys.argv) > 1 else '/data')/'x/data_for_paper'
OUT = Path(__file__).resolve().with_name('pool.json')
BOHR_A = 0.529177210903


def blocks(text):
    """Conductivity blocks in print order: ('SERTA', rows), ('iter1', rows), ..., final converged block last."""
    lines, out, label, i = text.split('\n'), [], None, 0
    while i < len(lines):
        line = lines[i]
        if 'self-energy relaxation time approximation' in line:
            label = 'SERTA'
        m = re.match(r'\s*Iteration number:\s+(\d+)', line)
        if m:
            label = 'iter' + m.group(1)
        if '[Ohm.cm]^-1' in line:
            j = i + 1
            while j < len(lines) and (not lines[j].strip() or lines[j].strip().startswith('=')):
                j += 1
            rows = []
            while j + 2 < len(lines) and len(lines[j].split()) == 7:
                T, ef, dos, pop, a, b, c = map(float, lines[j].split())
                r1, r2 = lines[j + 1].split(), lines[j + 2].split()
                sigma = [[a, b, c], [float(x) for x in r1[1:4]], [float(x) for x in r2[1:4]]]
                rows.append({'T_K': T, 'fermi_eV': ef, 'dos_states_per_Ry': dos, 'sigma': sigma})
                j += 3
            err = None
            k = j
            while k < len(lines) and not lines[k].strip():
                k += 1
            if k < len(lines) and 'Max error' in lines[k]:
                err = float(lines[k].split()[0])
            out.append((label, rows, err))
            label = None
            i = j
            continue
        i += 1
    return out


def ibte_params(text):
    p = {}
    for key in ('nkf1', 'nqf1', 'nk1', 'nq1', 'fsthick', 'temps', 'lifc', 'dvscf_dir', 'amass(1)', 'mp_mesh_k',
                'iterative_bte', 'nbndsub', 'degaussw'):
        m = re.search(r'^\s*%s\s*=\s*(.+?),?\s*$' % re.escape(key), text, re.M)
        if m:
            p[key] = m.group(1).strip().rstrip(',').strip("'")
    return p


def run_record(path, run_id):
    d = path.parent
    out = path.read_text(errors='ignore')
    tin = (d/'ibte.in').read_text()
    scf = (d/'scf.in').read_text()
    celldm = float(re.search(r'celldm\(1\)\s*=\s*([\d.]+)', scf).group(1))
    rec = {'id': run_id, 'metal': run_id.split('/')[0], 'archive_files': [str(path.relative_to(SRC.parent)),
           str((d/'ibte.in').relative_to(SRC.parent)), str((d/'scf.in').relative_to(SRC.parent))],
           'ibte_in': ibte_params(tin), 'lattice_parameter_bohr': celldm,
           'lattice_parameter_A': round(celldm * BOHR_A, 6)}
    km = re.findall(r'Using uniform MP k-mesh:\s+(\d+)\s+(\d+)\s+(\d+)', out)
    qm = re.findall(r'Using uniform q-mesh:\s+(\d+)\s+(\d+)\s+(\d+)', out)
    rec['fine_k_mesh'] = [int(x) for x in km[0]] if km else None
    rec['fine_q_mesh'] = [int(x) for x in qm[0]] if qm else None
    b = blocks(out)
    if not b:
        rec['results'] = None
        rec['note'] = 'archived output contains no transport results (the file repeats the input)'
        return rec
    serta = [x for x in b if x[0] == 'SERTA'][0]
    iters = [x for x in b if x[0].startswith('iter')]
    final = b[-1]                     # the last IBTE iteration, printed once the max error is below threshold
    if not final[0].startswith('iter') or final[2] is None or final[2] > 1e-5:
        raise ValueError(f'{run_id}: last block {final[0]} is not a converged IBTE iteration')
    results = []
    for n, row in enumerate(serta[1]):
        results.append({'T_K': row['T_K'], 'fermi_eV': row['fermi_eV'], 'dos_states_per_Ry': row['dos_states_per_Ry'],
                        'serta_sigma_per_ohm_cm': row['sigma'],
                        'ibte_iterations': [{'iteration': int(lab[4:]), 'sigma_xx_per_ohm_cm': rows[n]['sigma'][0][0],
                                             'max_error': err} for lab, rows, err in iters],
                        'ibte_final_sigma_per_ohm_cm': final[1][n]['sigma']})
        assert abs(final[1][n]['T_K'] - row['T_K']) < 1e-9
    rec['results'] = results
    rec['ibte_max_error_last'] = final[2]
    rec['n_ibte_iterations'] = len(iters)
    return rec


def ifc_record(path, ifc_id, kind, lattice_T, sampling_T=None):
    text = path.read_text()
    blk = lambda tag: re.search(r'<%s[^>]*>(.*?)</%s>' % (re.escape(tag), re.escape(tag)), text, re.S).group(1)
    nq = [int(x) for x in blk('MESH_NQ1_NQ2_NQ3').split()]
    C = [None] * (nq[0] * nq[1] * nq[2])
    for m in re.finditer(r'<s_s1_m1_m2_m3\.1\.1\.(\d+)\.(\d+)\.(\d+)>\s*<IFC>(.*?)</IFC>', text, re.S):
        i, j, k = (int(m.group(n)) - 1 for n in (1, 2, 3))
        C[(i * nq[1] + j) * nq[2] + k] = [float(x) for x in m.group(4).split()]
    assert all(c is not None and len(c) == 9 for c in C), ifc_id
    rec = {'id': ifc_id, 'metal': ifc_id.split(':')[0], 'kind': kind, 'lattice_temperature_K': lattice_T,
           'sampling_temperature_K': sampling_T, 'archive_file': str(path.relative_to(SRC.parent)),
           'header_alat_bohr': float(blk('CELL_DIMENSIONS').split()[0]),
           'header_volume_bohr3': float(blk('UNIT_CELL_VOLUME_AU')),
           'at_alat': [[float(x) for x in row.split()] for row in blk('AT').strip().split('\n')],
           'mass_amu': float(blk('MASS.1')), 'mesh': nq, 'C_Ry_per_bohr2': C}
    if kind == 'TDEP':
        uc = (path.parents[2]/'infile.ucposcar').read_text().split('\n')
        a1 = [float(x) for x in uc[2].split()]
        rec['tdep_cell_alat_bohr'] = round(2 * max(abs(x) for x in a1) / BOHR_A, 6)
        rec['tdep_iteration'] = path.parents[1].name
    return rec


def main():
    runs = []
    for path in sorted(SRC.glob('*/4-epw/**/ibte.out*')):
        rel = path.relative_to(SRC)
        if 'spec' in rel.parts or 'comp_EPW_interpolation' in rel.parts:
            continue
        parent = str(rel.parent)
        run_id = parent + '/' + path.name if path.name != 'ibte.out' else parent
        runs.append(run_record(path, run_id))
    ifcs = []
    temps = {'Pb': [0, 10, 110, 210, 310, 410, 510, 610], 'Al': [0, 481, 682, 900], 'Nb': [0, 1000, 1500, 2000]}
    for metal, Ts in temps.items():
        for T in Ts:
            ifcs.append(ifc_record(SRC/f'{metal}/2-ph/ph-{T}K/espresso.xml', f'{metal}:DFPT:a({T}K)', 'DFPT', T))
            if T:
                f = sorted(SRC.glob(f'{metal}/3-ph-anharm-TDEP/{T}K/sampling.*/iter.*/matdyn/updated_qe.fc.xml'))
                assert len(f) == 1, f
                ifcs.append(ifc_record(f[0], f'{metal}:TDEP{T}K:a({T}K)', 'TDEP', T, T))
    for T in (210, 410, 610):
        f = SRC/f'Pb/3-ph-anharm-TDEP/0K-alat/sampling.{T}k/iter.005/matdyn/updated_qe.fc.xml'
        ifcs.append(ifc_record(f, f'Pb:TDEP{T}K:a(0K)', 'TDEP', 0, T))
    readme = (SRC.parent/'data_for_paper/README.md') if (SRC.parent/'data_for_paper/README.md').exists() else None
    OUT.write_text(json.dumps({'source': 'Materials Cloud zxf0b-rs406 data_for_paper.tar.gz', 'runs': runs,
                               'ifcs': ifcs}) + '\n')
    print(len(runs), 'runs;', len(ifcs), 'force-constant sets ->', OUT, f'({OUT.stat().st_size / 1e6:.2f} MB)')
    for r in runs:
        print(r['id'], r['fine_k_mesh'], r['lattice_parameter_bohr'], 'no results' if r['results'] is None else
              [(x['T_K'], round(1e6 / x['ibte_final_sigma_per_ohm_cm'][0][0], 4)) for x in r['results']][:3])


if __name__ == '__main__':
    main()
