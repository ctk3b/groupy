import numpy as np
import pdb


class System():
    """System container class
    """
    def __init__(self):
        """
        """
        self.name = ''

        # numbas
        self.mol_id = int()

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

    def append_gbbs(self, gbbs):
        """
        """
        for gbb in gbbs:
            n_atoms = self.types.shape[0]
            assert n_atoms == self.xyz.shape[0]
            assert n_atoms == self.masses.shape[0]
            assert n_atoms == self.charges.shape[0]

            self.types = np.append(self.types, gbb.types)
            self.masses = np.append(self.masses, gbb.masses)
            self.charges = np.append(self.charges, gbb.charges)
            self.xyz = np.vstack((self.xyz, gbb.xyz))

            gbb.bonds[:, 1:] + n_atoms
            gbb.angles[:, 1:] + n_atoms
            gbb.dihedrals[:, 1:] + n_atoms
            self.bonds = np.vstack((self.bonds, gbb.bonds))
            self.angles = np.vstack((self.angles, gbb.angles))
            self.dihedrals = np.vstack((self.dihedrals, gbb.dihedrals))
