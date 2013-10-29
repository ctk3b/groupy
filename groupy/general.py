import numpy as np


# --- general handy stuff ---
def find_nearest(array, target):
    """Find index and value of array component whose numeric
    value is closest to 'target'
    """
    idx = np.abs(array - target).argmin()
    return idx, array[idx]


def anint(x):
    """
    """
    if x >= 0.5:
        return 1.
    elif x < -0.5:
        return -1.
    else:
        return 0.


def calc_angle(u, v, already_normed=True):
    """Compute the angle between 2 vectors
    """
    if already_normed:
        c = np.dot(u, v)
    else:
        c = np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)
    return np.arccos(c) * 180.0 / np.pi
