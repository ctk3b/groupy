import numpy as np

from groupy.monolayers import *
from groupy.mdio import *
from groupy.inertia import *

# --- general handy stuff ---
def find_nearest(array, target):
    """Find index and value of array component whose numeric
    value is closest to 'target'
    """
    idx = np.abs(array - target).argmin()
    return idx, array[idx]


def calc_com(xyz, mass=[]):
    """Calculate center of mass of a set of points
    """
    cum_mass = 0
    com = np.array([0., 0., 0.])
    if mass.size:
        assert len(xyz) == len(mass)
    else:
        mass = np.ones(xyz.shape)
    for i, coord in enumerate(xyz):
        com += coord * mass[i]
        cum_mass += mass[i]
    com /= cum_mass
    return com


