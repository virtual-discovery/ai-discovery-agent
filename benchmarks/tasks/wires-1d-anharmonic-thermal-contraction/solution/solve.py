"""Reference solution for wires-1d-anharmonic-thermal-contraction; runs as the Harbor oracle and in check_task."""
import json
import os
import re
import subprocess
import time
from pathlib import Path
import numpy as np

APP = Path(os.environ.get('APP_DIR', '/app'))
RY_EV = 13.605693122994
RY_CM = 109737.31570111268
KB_CM = 0.69503476
EV_J = 1.602176634e-19
BOHR_A = 0.529177210903
NUM = r'[-+]?\d+\.\d+(?:[eE][-+]?\d+)?'
evidence = []


def lab(endpoint, payload=None):
    cmd = ['python3', str(APP/'lab_client.py'), endpoint] + ([json.dumps(payload)] if payload is not None else [])
    out = subprocess.run(cmd, capture_output=True, text=True)
    result = json.loads(out.stdout)
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
    evidence.append(r['id'])
    return r['result']


# ---------------------------------------------------------------- dynamical matrices
def parse_dyn(text):
    lines = text.splitlines()
    h = lines[2].split()
    ntyp, nat, alat = int(h[0]), int(h[1]), float(h[3])
    cell = np.array([[float(x) for x in lines[4 + k].split()] for k in range(3)]) * alat
    i = 7
    mass = {}
    for k in range(ntyp):
        m = re.match(r"\s*(\d+)\s+'[^']*'\s+(\S+)", lines[i + k])
        mass[int(m.group(1))] = float(m.group(2))
    i += ntyp
    masses = np.array([mass[int(lines[i + k].split()[1])] for k in range(nat)])
    blocks = []
    while i < len(lines):
        if 'Matrix in cartesian axes' in lines[i]:
            j = i + 1
            while 'q = (' not in lines[j]:
                j += 1
            q = np.array([float(x) for x in lines[j].split('(')[1].split(')')[0].split()]) * 2 * np.pi / alat
            j += 1
            phi = np.zeros((nat, nat, 3, 3), complex)
            for _ in range(nat * nat):
                while not lines[j].strip():
                    j += 1
                a, b = (int(x) for x in lines[j].split())
                for r in range(3):
                    v = [float(x) for x in lines[j + 1 + r].split()]
                    phi[a - 1, b - 1, r] = np.array(v[0::2]) + 1j * np.array(v[1::2])
                j += 4
            blocks.append((q, phi))
            i = j
        else:
            i += 1
    return cell, masses, blocks


def dyn_set(files):
    names = sorted((n for n in files if 'Dynamical' in files[n]), key=lambda n: int(re.findall(r'(\d+)$', n)[0]))
    cell, masses, blocks = None, None, []
    for n in names:
        cell, masses, b = parse_dyn(files[n])
        blocks += b
    return cell, masses, blocks


def dmat(masses, phi):
    m = np.sqrt(np.repeat(masses, 3))
    n = len(masses)
    D = phi.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n) / np.outer(m, m)
    return 0.5 * (D + D.conj().T)


def freqs(masses, phi):
    w2 = np.linalg.eigvalsh(dmat(masses, phi))
    return np.sign(w2) * np.sqrt(np.abs(w2)) * RY_CM


def grid(files):
    """Per-q sorted frequencies (cm^-1) and the q vectors (1/bohr)."""
    cell, masses, blocks = dyn_set(files)
    return np.array([q for q, _ in blocks]), np.array([np.sort(freqs(masses, phi)) for _, phi in blocks])


def stiffness(files, ncell, fracs=(0.005, 0.01, 0.02)):
    """Long-wavelength LA stiffness M v^2 (eV per cell) from Fourier-interpolated real-space force constants."""
    cell, masses, blocks = dyn_set(files)
    c = cell[2, 2]
    ns = range(-(ncell // 2), ncell // 2 + 1)
    PhiR = {n: sum(phi * np.exp(1j * q[2] * n * c) for q, phi in blocks).real / len(blocks)
            * (0.5 if abs(n) == ncell // 2 else 1.0) for n in ns}
    ks = []
    for f in fracs:
        qz = f * 2 * np.pi / c
        D = dmat(masses, sum(P * np.exp(-1j * qz * n * c) for n, P in PhiR.items()))
        w2, v = np.linalg.eigh(D)
        zc = [np.sum(np.abs(v[2::3, k]) ** 2) for k in range(6)]
        k = int(np.argmax(zc))
        ks.append(masses.sum() * w2[k] / qz ** 2 * RY_EV)
    return float(np.mean(ks)), abs(np.linalg.det(cell)) * BOHR_A ** 3


def zero_modes(qs, W, n0):
    g = int(np.argmin(np.linalg.norm(qs, axis=1)))
    m = np.zeros(W.shape, bool)
    m[g, np.argsort(np.abs(W[g]))[:n0]] = True
    return m


def cv_ratio(F, T0, n0, deg=2):
    """C_V = T dS/dT of the SSCHA entropy: sum_modes c(w) (1 - T/w dw/dT), per vibrational mode in k_B."""
    temps = sorted(F)
    qs, W0 = F[T0]
    zm = zero_modes(qs, W0, n0)
    stack = np.array([F[T][1] for T in temps])
    tot = 0.0
    for iq in range(W0.shape[0]):
        for j in range(W0.shape[1]):
            if zm[iq, j]:
                continue
            w = W0[iq, j]
            x = w / (KB_CM * T0)
            c = x * x * np.exp(x) / np.expm1(x) ** 2
            dw = np.polyval(np.polyder(np.polyfit(temps, stack[:, iq, j], deg)), T0)
            tot += c * (1 - T0 / w * dw)
    return tot / (W0.size - n0)


def stress(excerpt, which):
    lines = excerpt.splitlines()
    for i, l in enumerate(lines):
        if which == 'sscha' and 'STRESS TENSOR' in l:
            return float(re.findall(NUM, lines[i + 3])[2])
        if which == 'average' and 'Ab initio average stress' in l:
            return float(re.findall(NUM, lines[i + 3])[2])
    raise ValueError('no stress')


def main():
    scf = buy({'kind': 'dft_file', 'wire': 'TaSe3', 'file': 'scf.out'})['text']
    p_line = re.findall(r'P=\s*(\S+)', scf)[-1]
    # TaSe3 harmonic stiffness from the 1D-ASR DFPT matrices (plain ph.x output violates the sum rules near Gamma)
    k, V = stiffness(buy({'kind': 'dfpt', 'wire': 'TaSe3', 'asr': '1d'})['files'], 6)
    # SSCHA stress along the wire (pressure-like sign convention, as in pw.x: the free-energy stress exceeds the
    # ensemble-average stress by a positive kinetic term that grows with T)
    logs = {T: buy({'kind': 'sscha_log', 'wire': 'TaSe3', 'temperature_K': T, 'configurations': 10000})
            for T in (0, 200, 300, 400)}
    P = {T: stress(r['final_population_excerpt'], 'sscha') for T, r in logs.items()}
    A = {T: stress(r['final_population_excerpt'], 'average') for T, r in logs.items()}
    kinetic = {T: P[T] - A[T] for T in P}
    Ts = [200, 300, 400]
    dPdT = np.polyfit(Ts, [P[T] for T in Ts], 1)[0]
    alpha = dPdT * 1e9 * V * 1e-30 / (k * EV_J) * 1e6
    # TaSe3 stability from the free-energy Hessian; the 4 rigid-body modes at Gamma (3 translations + axial
    # rotation, a 1D rotational sum-rule artefact) are not phonons
    soft = {}
    for T in (0, 500, 600, 800):
        qs, W = grid(buy({'kind': 'sscha_matrices', 'wire': 'TaSe3', 'temperature_K': T, 'configurations': 10000,
                          'matrix': 'hessian'})['files'])
        soft[T] = float(W[~zero_modes(qs, W, 4)].min())
    s = {T: np.sign(w) * w * w for T, w in soft.items()}
    T_inst = 600 + 200 * s[600] / (s[600] - s[800])
    # Heat capacity from the auxiliary (entropy-defining) frequencies of converged 10k runs bracketing 300 K
    cv, conv = {}, {}
    for wire, n0, extra_logs in (('TaSe3', 4, (100, 500, 600)), ('CuC2', 3, (100, 200, 300, 500, 600))):
        F = {T: grid(buy({'kind': 'sscha_matrices', 'wire': wire, 'temperature_K': T, 'configurations': 10000,
                          'matrix': 'auxiliary'})['files']) for T in (100, 200, 300, 500, 600)}
        conv[wire] = {T: 'satisfy the convergence criteria' in buy({'kind': 'sscha_log', 'wire': wire, 'temperature_K': T,
                                                                    'configurations': 10000})['final_population_excerpt']
                      for T in extra_logs}
        cv[wire] = cv_ratio(F, 300, n0)
    submission = {'TaSe3': {'axial_stiffness_eV': k, 'alpha_L_300K_1e-6_per_K': alpha, 'cv_ratio_300K': cv['TaSe3'],
                            'instability_temperature_K': T_inst},
                  'CuC2': {'cv_ratio_300K': cv['CuC2']}, 'evidence_ids': evidence}
    (APP/'submission.json').write_text(json.dumps(submission, indent=1) + '\n')
    (APP/'report.md').write_text(f"""# 1D wires: thermal contraction, heat capacity and stability (reference solution)

**Stress convention.** The SSCHA stress follows the pw.x convention (pressure-like, positive = compressive; the TaSe3
DFT reference prints P = {p_line} kbar). The free-energy stress exceeds the ab initio ensemble average by a positive
kinetic term growing with T ({', '.join(f'{T} K: {kinetic[T]:+.4f}' for T in sorted(kinetic))} GPa), which fixes the sign.
The average stress lacks this term and must not be differentiated.

**Thermal expansion.** At fixed length, alpha = V (dP_zz/dT)_L / k with V the DFT cell volume including vacuum
({V:.1f} A^3) and k = L0 d2E/dL2 per cell. dP_zz/dT (TaSe3, 200-400 K, converged 10k runs) = {dPdT:.3e} GPa/K.
k = {k:.1f} eV per cell from the long-wavelength LA branch of the 1D-ASR DFPT force constants; the plain ph.x
matrices violate the (1D) acoustic sum rules near Gamma, which blows up the small-q slope, and the single archived
strained calculation (+4.7 %) gives only an anharmonic finite-strain secant. alpha_L(300 K) = {alpha:.2f}e-6 /K:
TaSe3 contracts. The negative sign comes from the flexural and torsional modes (membrane/tension effect); the archive
contains only one strained calculation, so a quantitative QHA alpha cannot be reconstructed from it.

**Heat capacity.** The SSCHA entropy is the harmonic entropy of the auxiliary frequencies, so
C_V = sum c(w)(1 - T/w dw/dT); hardening frequencies lower C_V below the harmonic value, softening raises it.
Per vibrational mode (zero modes: 3 for the straight CuC2 chain, 4 for TaSe3 including the axial rotation):
CuC2 {cv['CuC2']:.3f}, TaSe3 {cv['TaSe3']:.3f} at 300 K (quadratic fits of each mode over 100-600 K; CuC2 400 K was
skipped because its minimization stopped at the Kong-Liu limit with residual centroid forces). Convergence flags of the
runs used: {conv}.

**Stability.** Auxiliary frequencies are positive by construction and cannot signal instabilities; the Hessian of the
positional free energy can. The small negative Hessian value at Gamma at all T is the axial rotation (1D rotational
invariance), not an instability. The softest genuine TaSe3 mode (q = 1/3 flexural pair) is {soft[600]:.1f} cm-1 at
600 K and {soft[800]:.1f} cm-1 at 800 K; interpolating omega^2 gives an onset near {T_inst:.0f} K (Hessians at 650 and
700 K are withheld). CuC2 Hessians show small imaginary, degeneracy-breaking flexural values that are within sampling
noise and are not treated as an instability. Stochastic errors (10k vs 100k configurations) limit the onset to
tens of kelvin.
""")
    lab('submit', submission)
    for _ in range(60):
        status = lab('submission')
        if status['status'] != 'validating':
            print(json.dumps(status))
            return
        time.sleep(2)


if __name__ == '__main__':
    main()
