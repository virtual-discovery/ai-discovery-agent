"""Numpy harmonic-phonon utilities for monatomic Bravais crystals from Quantum ESPRESSO q2r XML force constants.

Reads `espresso.xml` (q2r output) or TDEP's `updated_qe.fc.xml` (same schema), imposes the simple acoustic sum
rule, Fourier-interpolates with Wigner-Seitz image weights (as matdyn.x does) and returns frequencies in meV.
Only numpy is used, so the same code runs in the lab image (numpy 1.21).
"""
import re
import numpy as np

RY_MEV = 13605.693122994          # 1 Ry in meV
AMU_RY = 911.44424310865645       # 1 amu in Rydberg mass units (m_e / 2)
RY_CM1 = 109737.31568160          # 1 Ry in cm^-1
KB_MEV = 0.08617333262            # meV / K


def _block(text, tag):
    m = re.search(r'<%s[^>]*>(.*?)</%s>' % (re.escape(tag), re.escape(tag)), text, re.S)
    if not m:
        raise ValueError('tag %s not found' % tag)
    return m.group(1)


def read_ifc(path_or_text):
    text = path_or_text if '<Root>' in path_or_text[:200] else open(path_or_text).read()
    alat = float(_block(text, 'CELL_DIMENSIONS').split()[0])
    at = np.array(_block(text, 'AT').split(), float).reshape(3, 3)
    mass = float(_block(text, 'MASS.1'))
    nq = [int(x) for x in _block(text, 'MESH_NQ1_NQ2_NQ3').split()]
    C = np.zeros((nq[0], nq[1], nq[2], 3, 3))
    for m in re.finditer(r'<s_s1_m1_m2_m3\.1\.1\.(\d+)\.(\d+)\.(\d+)>\s*<IFC>(.*?)</IFC>', text, re.S):
        i, j, k = (int(m.group(n)) - 1 for n in (1, 2, 3))
        C[i, j, k] = np.array(m.group(4).split(), float).reshape(3, 3)
    found = len(re.findall(r'<s_s1_m1_m2_m3\.1\.1\.', text))
    if found != nq[0] * nq[1] * nq[2]:
        raise ValueError('expected %d IFC blocks, found %d' % (nq[0] * nq[1] * nq[2], found))
    vol = float(_block(text, 'UNIT_CELL_VOLUME_AU'))
    return {'alat_bohr': alat, 'at': at, 'mass_amu': mass, 'nq': nq, 'C': C, 'volume_bohr3': vol}


class Phonons:
    """Harmonic phonons of a monatomic crystal; q in Cartesian units of 2 pi / alat."""

    def __init__(self, ifc, asr=True):
        self.ifc = ifc
        at, nq, C = ifc['at'], ifc['nq'], ifc['C'].copy()
        if asr:   # simple ASR: on-site term = - sum of all others
            total = C.reshape(-1, 3, 3).sum(axis=0)
            C[0, 0, 0] -= total
        sup = np.array([nq[i] * at[i] for i in range(3)])            # supercell vectors (alat units)
        shifts = np.array([(a, b, c) for a in range(-2, 3) for b in range(-2, 3) for c in range(-2, 3)])
        R, W, F = [], [], []
        for i in range(nq[0]):
            for j in range(nq[1]):
                for k in range(nq[2]):
                    r0 = i * at[0] + j * at[1] + k * at[2]
                    imgs = r0 + shifts @ sup
                    d = np.linalg.norm(imgs, axis=1)
                    close = imgs[np.abs(d - d.min()) < 1e-6 * max(1.0, d.min())]
                    for r in close:
                        R.append(r)
                        W.append(1.0 / len(close))
                        F.append(C[i, j, k])
        self.R, self.W, self.F = np.array(R), np.array(W), np.array(F)
        self.mass = ifc['mass_amu'] * AMU_RY

    def dynmat(self, q):
        q = np.atleast_2d(q)
        phase = np.exp(-2j * np.pi * q @ self.R.T) * self.W               # (nq, nR)
        D = np.einsum('qr,rab->qab', phase, self.F) / self.mass
        return 0.5 * (D + np.conj(np.swapaxes(D, 1, 2)))

    def frequencies(self, q):
        """Frequencies in meV (negative for imaginary), sorted, shape (nq, 3)."""
        w2 = np.linalg.eigvalsh(self.dynmat(q))
        return np.sign(w2) * np.sqrt(np.abs(w2)) * RY_MEV

    def grid(self, n, chunk=4000):
        """Frequencies on a uniform n x n x n Gamma-centred grid of the full Brillouin zone."""
        f = np.arange(n) / n
        frac = np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T
        bg = np.linalg.inv(self.ifc['at']).T                               # reciprocal vectors, 2 pi / alat units
        qc = frac @ bg
        out = np.concatenate([self.frequencies(qc[s:s + chunk]) for s in range(0, len(qc), chunk)])
        return out


def moments(freqs_meV, cutoff_meV=0.5):
    """Phonon-DOS moments over all modes above a small cutoff (excludes the Gamma acoustic modes)."""
    w = freqs_meV.ravel()
    w = w[w > cutoff_meV]
    return {'mean_meV': float(w.mean()), 'inv1_meV-1': float(np.mean(1 / w)), 'inv2_meV-2': float(np.mean(1 / w ** 2)),
            'second_meV2': float(np.mean(w ** 2)), 'omega_log_meV': float(np.exp(np.mean(np.log(w)))),
            'max_meV': float(w.max()), 'n_imaginary_modes': int((freqs_meV.ravel() < -cutoff_meV).sum())}


def from_record(rec):
    """Build Phonons from a pool/lab force-constant record (flattened C on the q2r mesh)."""
    nq = rec['mesh']
    C = np.array(rec['C_Ry_per_bohr2'], float).reshape(nq[0], nq[1], nq[2], 3, 3)
    return Phonons({'alat_bohr': rec['header_alat_bohr'], 'at': np.array(rec['at_alat'], float),
                    'mass_amu': rec['mass_amu'], 'nq': nq, 'C': C, 'volume_bohr3': rec['header_volume_bohr3']})


def folded_grid(at, n):
    """Gamma-centred n^3 grid; returns q in Cartesian 2 pi/alat units folded into the first Brillouin zone."""
    bg = np.linalg.inv(np.array(at, float)).T
    f = np.arange(n) / n
    frac = np.array(np.meshgrid(f, f, f, indexing='ij')).reshape(3, -1).T
    qc = frac @ bg
    shifts = np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]) @ bg
    cand = qc[:, None, :] + shifts[None]
    return cand[np.arange(len(qc)), np.argmin(np.linalg.norm(cand, axis=2), axis=1)]


PATHS = {'fcc': [('G', (0, 0, 0)), ('X', (0, 1, 0)), ('W', (0.5, 1, 0)), ('L', (0.5, 0.5, 0.5)), ('G', (0, 0, 0)),
                 ('K', (0.75, 0.75, 0))],
         'bcc': [('G', (0, 0, 0)), ('H', (0, 0, 1)), ('N', (0.5, 0.5, 0)), ('G', (0, 0, 0)), ('P', (0.5, 0.5, 0.5)),
                 ('H', (0, 0, 1))]}


def lattice_type(at):
    at = np.abs(np.array(at, float))
    return 'fcc' if np.sum(at < 1e-8) == 3 else 'bcc'


def report(rec, n, eigenvectors=False, points_per_segment=40):
    """Phonon summary served by the lab for one force-constant record."""
    ph = from_record(rec)
    q = folded_grid(rec['at_alat'], n)
    if eigenvectors:
        w2, e = np.linalg.eigh(ph.dynmat(q))
        w = np.sign(w2) * np.sqrt(np.abs(w2)) * RY_MEV
    else:
        w = ph.frequencies(q)
    real = w[w > 0.05]
    nbin = 60
    top = float(np.ceil(w.max() * 1.02 * 10) / 10)
    hist, edges = np.histogram(real, bins=nbin, range=(0.0, top))
    mom = {'mean_meV': float(real.mean()), 'inverse_mean_meV-1': float(np.mean(1 / real)),
           'inverse_square_mean_meV-2': float(np.mean(real ** -2)), 'omega_log_meV': float(np.exp(np.mean(np.log(real)))),
           'rms_meV': float(np.sqrt(np.mean(real ** 2))), 'max_meV': float(w.max()),
           'modes_below_0.05_meV_excluded': int(w.size - real.size),
           'imaginary_modes': int((w < -0.05).sum())}
    kind = lattice_type(rec['at_alat'])
    labels, pts, path_q = [], PATHS[kind], []
    for (la, a), (lb, b) in zip(pts[:-1], pts[1:]):
        a, b = np.array(a, float), np.array(b, float)
        for t in np.linspace(0, 1, points_per_segment, endpoint=False):
            path_q.append(a + t * (b - a))
        labels.append(la)
    path_q.append(np.array(pts[-1][1], float))
    labels.append(pts[-1][0])
    path_q = np.array(path_q)
    out = {'lattice': kind, 'grid': {'n': n, 'description': 'Gamma-centred n x n x n grid of the full Brillouin zone; '
                                     'q in Cartesian units of 2 pi / alat, folded into the first Brillouin zone',
                                     'q_cart_2pi_over_alat': np.round(q, 6).tolist(),
                                     'frequencies_meV': np.round(w, 5).tolist()},
           'dos': {'bin_edges_meV': np.round(edges, 4).tolist(), 'mode_fraction': np.round(hist / real.size, 6).tolist()},
           'moments': mom,
           'dispersion': {'labels': labels, 'points_per_segment': points_per_segment,
                          'q_cart_2pi_over_alat': np.round(path_q, 5).tolist(),
                          'frequencies_meV': np.round(ph.frequencies(path_q), 4).tolist()}}
    if eigenvectors:
        out['grid']['eigenvectors_real_imag'] = [[[[round(float(c.real), 6), round(float(c.imag), 6)] for c in e[i][:, m]]
                                                  for m in range(3)] for i in range(len(q))]
        out['grid']['eigenvector_layout'] = '[q][mode][cartesian component][real, imag], modes ordered as frequencies'
    return out
