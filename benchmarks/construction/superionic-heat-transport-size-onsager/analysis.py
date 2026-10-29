"""Reference protocol and estimator variants on the served (0.08 ps, int16) series. Runs in the runtime image.

    docker run ... -v $PWD:/c -v <cache>:/cache ai-theorist-benchmark-runtime:20261007 python3 -B /c/analysis.py

Writes /cache/kappa_table.json: for every run, kappa (W/m/K) from
  mc      multicomponent cepstral (heat flux + anion velocity sum), full served band (6.25 THz), AIC   [reference]
  mc3     same, low-pass F* = 3 THz
  plain   single-component cepstral of the raw LAMMPS heat flux, full band
  plain3  same, F* = 3 THz
  direct  multicomponent direct Green-Kubo: L00 - L01^2/L11 from running integrals, averaged over 3-6 ps
with kappa = S_red(0) / (2 V kB <T>^2) for the extensive flux J*V of compute heat/flux, <T> the run's mean temperature.
"""
import json
import sys
import numpy as np
sys.path.insert(0, '/c')
import cepstral as C

DT = 0.08


def load(pool, meta, key):
    q = pool[key].astype(np.float64) * np.array(meta[key]['scale'])
    return q[:, :3], q[:, 3:]


def kappas(J, V, vol, T):
    pref = C.kappa_prefactor(1.0, T) / vol / 2.0
    f, S, nl = C.spectra([J, V], DT)
    k3 = int(round(3.0 / f[1]))
    red = C.reduced(S)
    out = {}
    for name, spec, m, cut in (('mc', red, 2, None), ('mc3', red, 2, k3), ('plain', S[:, 0, 0], 1, None),
                               ('plain3', S[:, 0, 0], 1, k3)):
        s0, e0, P, _ = C.cepstral(spec, 3, m, f_cut_index=cut)
        out[name] = s0 * pref
        out[name + '_err'] = e0 * pref
        out[name + '_P'] = P
    # direct multicomponent GK
    n = len(J)
    X = np.stack([J, V])                                   # (2, n, 3)
    Fx = np.fft.rfft(X, n=2 * n, axis=1)
    m = int(6.0 / DT) + 1
    acf = np.fft.irfft(np.einsum('ifl,jfl->ijf', np.conj(Fx), Fx), axis=2)[:, :, :m] / (n - np.arange(m)) / 3
    L = np.concatenate([np.zeros((2, 2, 1)), np.cumsum(0.5 * (acf[:, :, 1:] + acf[:, :, :-1]) * DT, axis=2)], axis=2)
    L01 = 0.5 * (L[0, 1] + L[1, 0])
    kt = (L[0, 0] - L01 ** 2 / L[1, 1]) * 2 * pref
    out['direct'] = float(kt[int(3.0 / DT):].mean())
    out['plain_direct'] = float((L[0, 0] * 2 * pref)[int(3.0 / DT):].mean())
    return out


def main():
    pool = np.load('/cache/pool.npz')
    meta = json.load(open('/cache/pool_meta.json'))
    table = {}
    for key in meta:
        J, V = load(pool, meta, key)
        r = kappas(J, V, meta[key]['volume_A3'], meta[key]['mean_temperature_K'])
        table[key] = {k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()}
        print(key, json.dumps(table[key]), flush=True)
    json.dump(table, open('/cache/kappa_table.json', 'w'), indent=0)


if __name__ == '__main__':
    main()
