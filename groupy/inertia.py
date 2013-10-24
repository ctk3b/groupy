import numpy as np


def calc_inertia_tensor(xyz, masses):
    """Calculates the moment of inertia tensor from a given set of points

    Args:
        xyz (np.ndarray):
        mass (np.ndarray):
    Returns:
        I (np.ndarray): moment of inertia tensor
    """
    I = np.zeros(shape=(3,3))

    com = calc_com(xyz, masses)
    for atom, coord0 in enumerate(xyz):
        mass = masses[atom]
        coord = coord0 - com
        I[0,0] += mass * (coord[1] * coord[1] + coord[2] * coord[2])
        I[1,1] += mass * (coord[0] * coord[0] + coord[2] * coord[2])
        I[2,2] += mass * (coord[0] * coord[0] + coord[1] * coord[1])
        I[0,1] -= mass * coord[0] * coord[1]
        I[0,2] -= mass * coord[0] * coord[2]
        I[1,2] -= mass * coord[1] * coord[2]
    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]
    return I

def calc_vector(I):
    """Calculates characteristic vector describing, e.g. a polymer

    Args:
        I (np.ndarray): 3x3 moment of inertia tensor
    Returns:
        direction (np.ndarray): characteristic vector
    """
    w, v = np.linalg.eig(I)
    direction = v[np.argmin(w)]
    return direction
