import numpy as np
from  scipy.spatial import cKDTree
import pdb
import math


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

def calc_distance_sq_pbc(point1, point2, box):
    """Naive squared distance calculation considering minimum image
    """
    assert point1.shape == point2.shape
    assert box.length.shape == point1.shape
    r_12 = 0.0
    for i, k in enumerate(point1):
        diff = k - point2[i]
        diff -= box.length[i] * anint(diff / box.length[i])
        r_12 += diff * diff
    return r_12

def calc_distance_pbc(x0, x1, dimensions):
    """Vectorized distance calculation considering minimum image
    """
    d = np.abs(x0 - x1)
    d = np.where(d > 0.5 * dimensions, dimensions - d, d)
    return np.sqrt((d ** 2).sum(axis=-1))
