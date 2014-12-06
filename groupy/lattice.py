import numpy as np

class Lattice():
    """Class to hold lattice points

    """

    def __init__(self):
        self.points = np.empty(shape=(0,3), dtype=float)


    def grid_mask_2d(self, n, m, center='origin', box=None):
        """ """
        mask = np.zeros(shape=(n*m, 3), dtype=float)
        for i in range(n):
            for j in range(m):
                mask[i*m + j, 0] = float(i) / float(n) + 0.5/n
                if center == 'origin':
                    mask[i*m + j, 0] -= 0.5
                mask[i*m + j, 1] = float(j) / float(m) + 0.5/m
                if center == 'origin':
                    mask[i*m + j, 1] -= 0.5


        # fit points nicely into box
        if box:
            mask = np.multiply(box.lengths, mask)
            mask[:] += box.mins - np.array(
                            [np.amin(mask[:, 0]),
                             np.amin(mask[:, 1]),
                             np.amin(mask[:, 2])]) + np.array(
                            [0.5*box.lengths[0]/float(n),
                             0.5*box.lengths[1]/float(m),
                             0.0])
            mask[:, 2] = 0.0
        self.points = mask

    def grid_mask_3d(self, n, m, l, center='origin', box=None):
        """ """
        mask = np.zeros(shape=(n*m*l, 3), dtype=float)
        for i in range(n):
            for j in range(m):
                for k in range(l):
                    mask[i*m*l + j*l + k, 0] = float(i) / n + 0.5/n
                    if center == 'origin':
                        mask[i*m*l + j*l + k, 0] -= 0.5
                    mask[i*m*l + j*l + k, 1] = float(j) / m + 0.5/m
                    if center == 'origin':
                        mask[i*m*l + j*l + k, 1] -= 0.5
                    mask[i*m*l + j*l + k, 2] = float(k) / l + 0.5/l
                    if center == 'origin':
                        mask[i*m*l + j*l + k, 2] -= 0.5
        import pdb
        if box:
            mask = np.multiply(box.lengths, mask)
            mask[:] += box.mins - np.array(
                            [np.amin(mask[:, 0]),
                             np.amin(mask[:, 1]),
                             np.amin(mask[:, 2])]) + np.array([
                             0.5*box.lengths[0]/float(n),
                             0.5*box.lengths[1]/float(m),
                             0.5*box.lengths[2]/float(l)])
        self.points = mask
