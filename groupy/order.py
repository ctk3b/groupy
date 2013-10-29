import numpy as np
import pdb

from groupy.general import *


def calc_director(I):
    """Calculates characteristic vector describing, e.g. a polymer

    Args:
        I (np.ndarray): 3x3 moment of inertia tensor
    Returns:
        director (np.ndarray): characteristic vector
    """
    w, v = np.linalg.eig(I)
    director = v[:, np.argmin(w)]
    return director


def calc_Q_tensor(directors):
    """
    """
    normed = directors / np.sqrt((directors ** 2).sum(-1))[..., np.newaxis]
    Q = np.zeros(shape=(3,3))
    for vector in normed:
        Q[0,0] += 3.0 * vector[0] * vector[0] - 1
        Q[0,1] += 3.0 * vector[0] * vector[1]
        Q[0,2] += 3.0 * vector[0] * vector[2]
        Q[1,0] += 3.0 * vector[1] * vector[0]
        Q[1,1] += 3.0 * vector[1] * vector[1] - 1
        Q[1,2] += 3.0 * vector[1] * vector[2]
        Q[2,0] += 3.0 * vector[2] * vector[0]
        Q[2,1] += 3.0 * vector[2] * vector[1]
        Q[2,2] += 3.0 * vector[2] * vector[2] - 1
    Q /= (2.0 * normed.shape[0])
    return Q


def calc_S2(Q):
    """Calculate nematic order parameter (S2)
    """
    w, v = np.linalg.eig(Q)
    S2 = w.max()
    director = v[:, np.argmin(w)]
    return S2, director
