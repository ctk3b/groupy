import numpy as np


# --- general handy stuff ---
def find_nearest(array, target):
    """Find index and value of array component whose numeric
    value is closest to 'target'
    """
    idx = np.abs(array - target).argmin()
    return idx, array[idx]


def calc_com(xyz, mass):
    """Calculate center of mass of a set of points
    """
    assert xyz.shape[0] == mass.shape[0]

    cum_mass = 0
    com = np.array([0., 0., 0.])
    for i, coord in enumerate(xyz):
        com += coord * mass[i]
        cum_mass += mass[i]
    com /= cum_mass
    return com


