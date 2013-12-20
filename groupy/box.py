import numpy as np
import pdb

class Box():
    """Class to hold box information.

    TODO:
        -functions to update all box properties when one is changed
    """
    def __init__(self, t1, t2=None, t3=None):
        """
        """
        self.mins = np.zeros(3, dtype='f')
        self.maxs = np.zeros(3, dtype = 'f')
        self.length = np.zeros(3, dtype = 'f')
        # blank constructor
        if t1 is None:
            pass
        # box lengths passed as arguments
        elif type(t1) is int or type(t1) is float:
            if t2 is None:
                t2 = t1
            if t3 is None:
                t3 = t1
            vals = np.array([t1, t2, t3])
            self.maxs = vals
            self.length = vals
        # 3D box dimensions passed as arguments
        elif t1.shape[0] == 3:
            self.maxs = t1
            if t2 is None:
                pass
            else:
                self.mins = t2
            if t3 is None:
                self.length = t1 - t2
        self.dims = np.vstack([self.mins, self.maxs]).T
        self.volume = self.length[0] * self.length[1] * self.length[2]

    def update(self, changed):
        """
        """
        if changed == 'dims':
            self.mins = self.dims[:, 0]
            self.maxs = self.dims[:, 1]
            self.length = self.maxs - self.mins
            self.volume = self.length[0] * self.length[1] * self.length[2]
