import numpy as np
from scipy.spatial import cKDTree

from box import *

import pdb


class System():
    """System container class
    """
    def __init__(self, box=Box(0)):
        """
        """
        self.name = ''

        # numbas
        self.resids = np.empty(shape=0, dtype='int')
        self.resnames = np.empty(shape=0, dtype='str')

        # per atom
        self.types = np.empty(shape=(0), dtype='int')
        self.masses = np.empty(shape=(0))
        self.charges = np.empty(shape=(0))
        self.xyz = np.empty(shape=(0, 3))

        # connectivity
        self.bonds = np.empty(shape=(0, 3), dtype='int')
        self.angles = np.empty(shape=(0, 4), dtype='int')
        self.dihedrals = np.empty(shape=(0, 5), dtype='int')
        self.impropers = np.empty(shape=(0, 5), dtype='int')

        # forcefield info
        self.pair_types = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.dihedral_types = dict()
    
        self.box = box 

    def append_gbbs(self, gbbs):
        """
        """
        for gbb in gbbs:
            n_new_atoms = gbb.types.shape[0]
            assert n_new_atoms == gbb.xyz.shape[0]
            if gbb.masses.shape[0] > 0:
                assert n_new_atoms == gbb.masses.shape[0]
            if gbb.charges.shape[0] > 0:
                assert n_new_atoms == gbb.charges.shape[0]

            n_sys_atoms = self.types.shape[0]

            self.types = np.append(self.types, gbb.types)
            self.masses = np.append(self.masses, gbb.masses)
            self.charges = np.append(self.charges, gbb.charges)
            self.xyz = np.vstack((self.xyz, gbb.xyz))

            gbb.bonds[:, 1:] + n_sys_atoms
            gbb.angles[:, 1:] + n_sys_atoms
            gbb.dihedrals[:, 1:] + n_sys_atoms
            self.bonds = np.vstack((self.bonds, gbb.bonds))
            self.angles = np.vstack((self.angles, gbb.angles))
            self.dihedrals = np.vstack((self.dihedrals, gbb.dihedrals))

    def init_atom_kdtree(self):
        """
        """
        self.atom_kdtree = cKDTree(self.xyz)

    def get_atoms_in_range(self, point, radius, max_items=50):
        """
        """
        if not hasattr(self, 'atom_kdtree'):
            self.init_atom_kdtree()

        distances, indices = self.atom_kdtree.query(point, max_items)

        neighbors = []
        for index, distances in zip(indices, distances):
            if distance <= radius:
                neighbors.append(index)
            else:
                break

        return neighbors
