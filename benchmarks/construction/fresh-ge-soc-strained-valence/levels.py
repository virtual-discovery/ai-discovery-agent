"""Gamma-point level analysis shared by the construction scripts (same estimators as solution/solve.py).

Spin-orbit runs: 36 spinor bands; at Gamma bands 21-22 are the s-like bonding state (Gamma6+), 23-30 hold the
s-like antibonding Gamma7- (inverted below the p manifold in PBE) and the p-like Gamma7+/Gamma8+ manifold, in an
order that changes with strain. Identification is by the bands.x parity labels when present, otherwise by pairing
levels to the labelled reference run.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows(path=None):
    out = {}
    for name in ([path] if path else ['exploration_runs.jsonl', 'reference_runs.jsonl']):
        p = HERE/name
        if p.exists():
            for line in p.read_text().splitlines():
                if line.strip():
                    r = json.loads(line)
                    out[r['tag']] = r
    return out


def pairs(e, lo=22, hi=30):
    """Kramers-pair energies (eV) of spinor bands lo+1..hi at one k point (1-based band indices)."""
    return [round((e[i] + e[i + 1]) / 2, 5) for i in range(lo, hi, 2)]


def gamma_levels_soc(result):
    """{'s': E(odd-parity pair), 'p': [three even-parity pair energies, high to low]} at Gamma from a labelled run."""
    band = next(b for b in result['band_kpoints'] if b['k_cartesian_inverse_A'] == [0.0, 0.0, 0.0])
    s, p = None, []
    for lev in band['symmetry']['levels']:
        i, j = lev['bands']
        if 23 <= i <= 30:
            n = j - i + 1
            if lev['label'].split()[0].endswith('-'):
                s = lev['energy_eV']
            else:
                p += [lev['energy_eV']] * (n // 2)
    return {'s': s, 'p': sorted(p, reverse=True), 'labels': [(l['bands'], l['energy_eV'], l['label'])
                                                             for l in band['symmetry']['levels'] if 21 <= l['bands'][0] <= 32]}
