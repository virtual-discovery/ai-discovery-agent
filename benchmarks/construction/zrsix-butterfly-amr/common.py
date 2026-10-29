"""Shared settings of the ZrSiX construction (energy windows, Btau grid, field-direction menu)."""
import numpy as np

COMPOUNDS = ['ZrSiS', 'ZrSiSe', 'ZrSiTe']
GROUP = {'ZrSiS': [11, 12], 'ZrSiSe': [11, 12], 'ZrSiTe': [11, 12, 13]}       # 0-based band indices
# transport window of start energies (eV, model scale): [lo, hi]; chemical potentials served are lo+0.02 .. hi-0.02
WINDOW = {'ZrSiS': (9.668, 9.745), 'ZrSiSe': (8.882, 8.945), 'ZrSiTe': (8.830, 8.920)}
EBIN = 0.001
XGRID = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20.0])
MESH_TOP = 160
SUBMESH = [40, 80, 160]


def directions(name):
    """Field-direction menu: list of (plane, angle_deg)."""
    out = []
    if name == 'ZrSiSe':
        out += [('bc', float(t)) for t in range(0, 91, 5)]
        out += [('ac', float(t)) for t in (30, 45, 60, 90)]
        out += [('ab', 45.0)]
    elif name == 'ZrSiTe':
        out += [('bc', float(t)) for t in (-90, -80, -70, -60, -50, 0, 30, 60, 90)]
    else:
        out += [('bc', float(t)) for t in (0, 20, 40, 50, 60, 90)]
        out += [('ab', 45.0)]
    return out


def bvec(plane, ang):
    t = np.radians(ang)
    if plane == 'bc':
        return np.array([0.0, np.sin(t), np.cos(t)])
    if plane == 'ac':
        return np.array([np.sin(t), 0.0, np.cos(t)])
    if plane == 'ab':
        return np.array([np.cos(t), np.sin(t), 0.0])
    raise ValueError(plane)


def dkey(plane, ang):
    return f'{plane}_{ang:+06.1f}'

# electron-count service: chemical potentials allowed in COUNT_RANGE; histograms extend 0.3 eV beyond it
COUNT_RANGE = {'ZrSiS': (8.70, 10.45), 'ZrSiSe': (8.05, 9.75), 'ZrSiTe': (7.95, 9.70)}
COUNT_BIN = 0.0002
COUNT_MESHES = [40, 80, 120, 160]
# chemical potentials accepted by electron_count (inside the histogram span, which is COUNT_RANGE +- 0.3 eV)
COUNT_SERVE = {'ZrSiS': (8.70, 10.70), 'ZrSiSe': (8.05, 10.00), 'ZrSiTe': (7.95, 9.95)}
