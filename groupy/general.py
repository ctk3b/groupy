import numpy as np
from  scipy.spatial import cKDTree
import pdb
import math
import copy


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

def can_add(gbb, gbb_list, box, r_cut):
    """
    """
    if len(gbb_list) > 50:
        coms = np.asarray([i.com for i in gbb_list])
        d_coms = calc_distance_pbc(coms, gbb.com, box.maxs)
        test_list = [gbb_list[i] for i, d in enumerate(d_coms) if d < 20]
    else:
        test_list = gbb_list
    for test_atom in gbb.xyz:
        for item in test_list:
            r_ij = calc_distance_pbc(item.xyz, test_atom, box.maxs)
            if (r_ij < r_cut).any():
                return False
    return True

def add_to_box(gbb, gbb_list, n, box, dims=[True, True, True], r_cut=2.0, name='gbb'):
    """
    """
    added = list()
    while (len(added) < n):
        t_gbb = copy.deepcopy(gbb)
        t_gbb.shift_com_to_origin()
        t_gbb.rotate(180 * np.random.rand(3))
        x = np.random.uniform(box.mins[0], box.maxs[0])
        y = np.random.uniform(box.mins[1], box.maxs[1])
        z = np.random.uniform(box.mins[2], box.maxs[2])
        t_gbb.translate([x, y, z])
        t_gbb.calc_com()
        if can_add(t_gbb, added + gbb_list, box, r_cut):
            t_gbb.wrap(box, dims)
            added.append(t_gbb)
            if len(added) % 100 == 0:
                print "Added {0} #{1}".format(name, len(added))
    return added

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
