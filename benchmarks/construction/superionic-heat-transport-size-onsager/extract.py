"""Reduce the fetched heat-transport logs to the served format (run in the runtime image).

    docker run --rm --cpuset-cpus 6,7,22,23 --user $(id -u):$(id -g) -v $PWD:/c -v ~/.cache/ai-theorist/<id>:/cache \
        ai-theorist-benchmark-runtime:20261007 python3 -B /c/extract.py

Input:  /cache/raw/<material>/<set>/<n>x<n>x<n>/T_<T>_K/{nve,nvt}_log.npy  (from fetch.py; 100001 rows, 8 fs)
        plus each run's LAMMPS input deck.
Output: /cache/pool.npz      int16 arrays (n_frames, 6) per run: J_x, J_y, J_z, V_x, V_y, V_z
        /cache/pool_meta.json per-run metadata (scales, cell, volume, mean/std temperature, energy drift, deck)

Transformation: the first 100000 frames (800 ps) are band-limited exactly to the Nyquist frequency of a
0.08 ps grid (FFT truncation, i.e. an ideal low-pass at 6.25 THz) and resampled to 10000 frames; each column is
then stored as int16 with a per-column scale (quantisation noise < 1e-8 of the signal variance).
"""
import glob
import json
import numpy as np

RAW = '/cache/raw'
SETS = {('UO2', 'NVE'): 'UO2_heat/NVE_L_T', ('UO2', 'NVT'): 'UO2_heat/NVT_L_T_heat',
        ('CaF2', 'NVE'): 'CaF2/NVE_L_T_heat', ('CaF2', 'NVT'): 'CaF2/NVT_L_T_heat',
        ('PbF2', 'NVE'): 'PbF2/NVE_L_T_6p056_heat', ('PbF2', 'NVT'): 'PbF2/NVT_L_T_6p056_heat'}
ATOMS_PER_CELL = 12
DT_IN, DT_OUT, N_IN = 0.008, 0.08, 100000
N_OUT = N_IN * DT_IN / DT_OUT


def decimate(x):
    n_out = int(round(N_OUT))
    F = np.fft.rfft(x[:N_IN], axis=0)
    keep = n_out // 2 + 1
    G = F[:keep].copy()
    G[-1] *= 0.5 if n_out % 2 == 0 else 1.0       # split the Nyquist bin (real signal)
    return np.fft.irfft(G, n=n_out, axis=0) * (n_out / N_IN)


def main():
    arrays, meta = {}, {}
    for (mat, ens), folder in SETS.items():
        for f in sorted(glob.glob(f'{RAW}/{folder}/*/T_*_K/*_log.npy')):
            parts = f.split('/')
            n = int(parts[-3].split('x')[0])
            T = int(parts[-2].split('_')[1])
            a = np.load(f).astype(np.float64)
            assert len(a) >= N_IN + 1, (f, len(a))
            assert np.allclose(np.diff(a[:200, 0]), DT_IN, atol=1e-4), f
            x = decimate(a[:, 6:12])
            scale = np.abs(x).max(0) / 32000.0
            q = np.round(x / scale).astype(np.int16)
            key = f'{mat}-{ens}-{n}x{n}x{n}-{T}K'
            arrays[key] = q
            deck = open(glob.glob(f'{RAW}/{folder}/{parts[-3]}/{parts[-2]}/*.in')[0]).read()
            t = a[:N_IN, 0] - a[0, 0]
            meta[key] = {'material': mat, 'ensemble': ens, 'supercell': n, 'temperature_K': T,
                         'atoms': ATOMS_PER_CELL * n ** 3, 'Lx_A': float(a[0, 5]),
                         'volume_A3': float(a[0, 5]) ** 3,
                         'mean_temperature_K': float(a[:N_IN, 3].mean()),
                         'std_temperature_K': float(a[:N_IN, 3].std()),
                         'mean_pressure_bar': float(a[:N_IN, 4].mean()),
                         'total_energy_drift_eV_per_ns': float(np.polyfit(t, a[:N_IN, 1], 1)[0] * 1000),
                         'mean_potential_energy_eV': float(a[:N_IN, 2].mean()),
                         'scale': scale.tolist(), 'deck': deck}
            print(key, flush=True)
    np.savez('/cache/pool.npz', **arrays)
    json.dump(meta, open('/cache/pool_meta.json', 'w'), indent=0)
    print(len(arrays), 'runs')


if __name__ == '__main__':
    main()
