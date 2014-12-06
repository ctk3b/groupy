import numpy as np
import pdb

class Box():
    """Class to hold box information.

    TODO:
        -functions to update all box properties when one is changed
    """
    def __init__(self, lengths=None, mins=None, maxs=None, center=None):
        """Initialize a box.

        """
        if np.any(lengths):
            self.lengths = lengths
            self.mins = np.array(
                    [-lengths[0]/ 2,
                    -lengths[1]/2, 
                    -lengths[2]/2])
            self.maxs = np.array(
                    [lengths[0]/ 2,
                     lengths[1]/2,
                     lengths[2]/2])
                    
        # TODO: make this better
        elif np.any(mins) and np.any(maxs):
            self.mins = mins
            self.maxs = maxs
            self.lengths = np.array(
                   [maxs[0] - mins[0],
                    maxs[1] - mins[1],
                    maxs[2] - mins[2]])

    def update(self, changed):
        """
        """
        if changed == 'dims':
            self.mins = self.dims[:, 0]
            self.maxs = self.dims[:, 1]
            self.length = self.maxs - self.mins
            self.volume = self.length[0] * self.length[1] * self.length[2]
