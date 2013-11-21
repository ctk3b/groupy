import numpy as np

class Box():
    """Class to hold box information.
        
    """
    
    def __init__(self, t1, t2=0, t3=0):
        self.mins = np.zeros(3, dtype='f')
        self.maxs = np.zeros(3, dtype = 'f')
        self.length = np.zeros(3, dtype = 'f')
        if not t2:
            t2 = t1
        if not t3:
            t3 = t1
        vals = np.array([t1, t2, t3])
        self.maxs = vals
        self.length = vals
        self.dims = np.vstack([self.mins, self.maxs]).T

