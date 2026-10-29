"""Readers for the Grillo et al. archive and the estimators shared by make_task.py and write_design.py.

All energies are eV. Band indices of the QP tables follow the archive (Yambo/QE numbering of the
spinor calculation); k-point n of the QP tables sits at k = (n - 1)/48 in units of 2*pi/c.
"""
import re
from pathlib import Path
import numpy as np
from scipy.special import k0

CACHE = Path.home()/'.cache/ai-theorist/gwbse-1d-chain-excitons/materials_cloud'
CHAINS = ['S3', 'Te3', 'As2S3', 'Bi2Te3']
PREFIX = {'S3': 's3', 'Te3': 'te3', 'As2S3': 'as2s3', 'Bi2Te3': 'bi2te3'}
SC_FILE = {'S3': 'G3W3', 'Te3': 'G3W3', 'As2S3': 'G4W4', 'Bi2Te3': 'G3W3'}
VBM = {'S3': 18, 'Te3': 48, 'As2S3': 48, 'Bi2Te3': 78}            # QP-table numbering (= SOC band file)
VBM_NOSOC = {'S3': 9, 'Te3': 36, 'As2S3': 24, 'Bi2Te3': 51}        # numbering of the no-SOC band file
NK = 25
BRIGHT_FRACTION = 0.10
HA, BOHR = 27.211386, 0.529177
P_ENVIRON = 0.6798e-08   # Ry/bohr^3, printed by Environ for env_pressure = 0.0001 GPa


def read_header(path):
    """Comment header of a Yambo output without the ASCII-art banner (kept from the 'Version' line on)."""
    lines = [l.rstrip() for l in Path(path).read_text().splitlines() if l.startswith('#')]
    start = next(i for i, l in enumerate(lines) if 'Version' in l)
    return lines[start:]


def read_qp(chain, level):
    name = 'G0W0' if level == 'G0W0' else SC_FILE[chain]
    path = CACHE/f'GW/{chain}/o-{name}.qp'
    data = np.loadtxt(path, comments='#')
    return {'file': f'GW/{chain}/o-{name}.qp', 'header': read_header(path), 'rows': data}


def read_bse(chain, run):
    path = CACHE/f'GW/{chain}/o-{run}.alpha_q1_haydock_bse'
    header = read_header(path)
    cols = next(l for l in header if 'E[1]' in l)
    names = re.split(r'\s{2,}', cols.lstrip('#').strip())
    return {'file': f'GW/{chain}/o-{run}.alpha_q1_haydock_bse', 'header': header, 'columns': names,
            'data': np.loadtxt(path, comments='#')}


def read_bands(chain, soc):
    path = CACHE/f"GW/{chain}/{PREFIX[chain]}_{'soc' if soc else 'nosoc'}.bands.dat.gnu"
    blocks, cur = [], []
    for line in path.read_text().splitlines():
        s = line.split()
        if not s:
            if cur:
                blocks.append(cur)
                cur = []
            continue
        cur.append([float(x) for x in s])
    if cur:
        blocks.append(cur)
    a = np.array(blocks)
    k = a[0, :, 0]
    return {'file': str(path.relative_to(CACHE)), 'k_frac': 0.5 * k / k[-1], 'energies': a[:, :, 1]}


def read_screening(chain):
    eps_text = (CACHE/f'IR/{chain}/ph.epsil.out').read_text()
    block = eps_text[eps_text.rfind('Dielectric constant in cartesian axis'):].splitlines()[2:5]
    eps = [[float(x) for x in re.findall(r'-?\d+\.\d+', line)] for line in block]
    def cell(path):
        lines = Path(path).read_text().splitlines()
        i = next(j for j, l in enumerate(lines) if l.strip().upper().startswith('CELL_PARAMETERS'))
        return [[float(x) for x in lines[i + r].split()[:3]] for r in (1, 2, 3)]
    env = (CACHE/f'volume/{chain}/scf.environ.out').read_text()
    pv = float(re.findall(r'PV energy\s*=\s*([-\d.Ee+]+)\s*Ry', env)[-1])
    pressure = float(re.findall(r'external pressure in inter\. units\s*=\s*([-\d.Ee+]+)', env)[-1])
    nat = int(re.findall(r'nat\s*=\s*(\d+)', (CACHE/f'IR/{chain}/scf.in').read_text())[0])
    return {'dfpt_cell_A': cell(CACHE/f'IR/{chain}/scf.in'), 'epsilon_infinity_supercell': eps,
            'environ_cell_A': cell(CACHE/f'volume/{chain}/scf.in'), 'environ_pv_energy_Ry': pv,
            'environ_pressure_Ry_per_bohr3': pressure, 'atoms_per_cell': nat}


# ---------------------------------------------------------------- estimators
def band(rows, b, qp=True):
    k, bb = rows[:, 0].astype(int), rows[:, 1].astype(int)
    e = rows[:, 2] + (rows[:, 3] if qp else 0)
    return np.array([e[(k == i) & (bb == b)][0] for i in range(1, NK + 1)])


def qp_direct_gap(chain, level):
    rows = read_qp(chain, level)['rows']
    g = band(rows, VBM[chain] + 1) - band(rows, VBM[chain])
    return float(g.min()), int(g.argmin()) + 1


def local_maxima(E, y):
    return [i for i in range(1, len(E) - 1) if y[i] > y[i - 1] and y[i] >= y[i + 1]]


def first_bright(chain, run, fraction=BRIGHT_FRACTION, column=1):
    d = read_bse(chain, run)['data']
    E, y = d[:, 0], d[:, column]
    peaks = [i for i in local_maxima(E, y) if y[i] >= fraction * y.max()]
    return float(E[peaks[0]])


def ip_first_peak(chain, run, fraction=0.01):
    return first_bright(chain, run, fraction, column=3)


def ip_onset(chain, run, fraction=0.01):
    d = read_bse(chain, run)['data']
    return float(d[np.argmax(d[:, 3] > fraction * d[:, 3].max()), 0])


def dft_min_direct_gap(chain, soc):
    b = read_bands(chain, soc)
    v = VBM[chain] if soc else VBM_NOSOC[chain]
    g = b['energies'][v] - b['energies'][v - 1]
    return float(g.min())


def soc_gap_change(chain):
    """SOC change of the minimum direct DFT gap (SOC minus no SOC; each structure's own minimum)."""
    return dft_min_direct_gap(chain, True) - dft_min_direct_gap(chain, False)


def effective_masses(chain, level, kidx):
    """Electron and hole masses (m_e) from a 5-point parabola at QP-table k-point kidx (mirror at 0 and pi/c)."""
    rows = read_qp(chain, level)['rows']
    c_bohr = CHAIN_C_A[chain] / BOHR
    dk = 2 * np.pi / (48 * c_bohr)
    out = []
    for b, sign in ((VBM[chain] + 1, 1), (VBM[chain], -1)):
        E = band(rows, b)
        idx = [min(abs(j), 2 * (NK - 1) - abs(j)) for j in range(kidx - 3, kidx + 2)]
        p = np.polyfit(np.arange(-2, 3) * dk, E[idx] / HA, 2)
        out.append(sign / (2 * p[0]))
    return out


def wire_polarizability(chain):
    """Polarizability per unit length (bohr^2) from the DFPT supercell epsilon_zz: (eps_zz - 1) A / 4 pi."""
    s = read_screening(chain)
    a = np.array(s['dfpt_cell_A'])
    area = abs(np.linalg.det(a[:2, :2])) / BOHR ** 2
    return (s['epsilon_infinity_supercell'][2][2] - 1) * area / (4 * np.pi)


def quantum_radius(chain):
    s = read_screening(chain)
    volume = s['environ_pv_energy_Ry'] / s['environ_pressure_Ry_per_bohr3']
    return float(np.sqrt(volume / (CHAIN_C_A[chain] / BOHR) / np.pi))


def wire_exciton_binding(mu, alpha, r0, L=3000.0, N=1201):
    """Lowest bound state of the 1D screened-wire exciton W(q) = v/(1 + alpha q^2 v), v = 2 K0(|q| r0) (a.u.)."""
    k = (np.arange(N) - N // 2) * 2 * np.pi / L
    q = np.abs(k[:, None] - k[None, :])
    q[q == 0] = 1e-12
    v = 2 * k0(q * r0)
    W = v / (1 + alpha * q ** 2 * v)
    h = 2 * np.pi / L
    W[np.arange(N), np.arange(N)] = 2 * (-np.log(h * r0 / 4) - 0.5772156649 + 1)   # cell average of the log singularity
    H = np.diag(k ** 2 / (2 * mu)) - W * h / (2 * np.pi)
    return float(-np.linalg.eigvalsh(H)[0] * HA)


CHAIN_C_A = {}


def _lattice_constants():
    for chain in CHAINS:
        CHAIN_C_A[chain] = read_screening(chain)['dfpt_cell_A'][2][2] if CACHE.exists() else None


_lattice_constants()
