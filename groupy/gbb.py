import numpy as np
import copy

from groupy.mdio import *
from groupy.general import *


class gbb():
    """
    """
    def __init__(self):
        """
        """
        self.n_atoms = None

        self.types = None
        self.elements = None
        self.mol_id = None

        self.xyz = None
        self.bonds = None
        self.angles = None
        self.dihedrals = None

        self.masses = None
        self.charges = None
        self.com = None


    # --- calculable properties ---
    def calc_com(self):
        """Calculate center of mass of a set of points
        """
        assert self.xyz.shape[0] == self.masses.shape[0]
        cum_mass = 0
        com = np.array([0., 0., 0.])
        for i, coord in enumerate(self.xyz):
            com += coord * self.masses[i]
            cum_mass += self.masses[i]
        com /= cum_mass
        self.com = com


    def calc_inertia_tensor(self):
        """Calculates the moment of inertia tensor of the gbb

        Returns:
            I (np.ndarray): moment of inertia tensor
        """
        assert self.xyz.shape[0] == self.masses.shape[0]
        I = np.zeros(shape=(3,3))
        self.calc_com()
        for i, coord0 in enumerate(self.xyz):
            mass = self.masses[i]
            coord = coord0 - self.com
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

    
    def calc_r_gyr():
        pass

    # --- deformations ---
    def translate():
        pass
    def rotate():
        pass
    def scale():
        pass

    def unwrap(self, box, dim=[True, True, True]):
        """Unwrap periodic boundary conditons of an object

        Requires that object being unwrapped does not span more than half the 
        box length.
        """
        for i, atom in enumerate(self.xyz):
            if i != 0:
                for k in range(3):  # TODO: vectorize
                    if dim[k] == True:
                        dr = self.xyz[i, k] - self.xyz[0, k]
                        box_length = np.diff(box[k])
                        self.xyz[i, k] -= box_length * anint(dr / box_length)

    def wrap():
        pass

    # --- io ---
    def load_mass(self, file_name):
        self.masses = np.loadtxt(file_name) 


    def load_xyz(self, file_name):
        self.xyz, self.types = read_xyz(file_name) 
        self.n_atoms = self.xyz.shape[0]

