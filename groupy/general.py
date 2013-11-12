import numpy as np
from  scipy.spatial import cKDTree
import pdb


# --- general handy stuff ---
def find_nearest(array, target):
    """Find array component whose numeric value is closest to 'target'.
    """
    idx = np.abs(array - target).argmin()
    return idx, array[idx]


def anint(x):
    """Switching function used for PBC unwrapping.
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


# spatial search 
def get_points_in_range(array, point, radius, max_items=50):
    """
    """
    tree = cKDTree(array)

    distances, indices = tree.query(point, max_items)

    neighbors = []
    for index, distance in zip(indices, distances):
        if distance <= radius:
            neighbors.append(index)
        else:
            break

    return neighbors
