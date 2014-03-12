import numpy as np
from scipy.spatial import cKDTree

from box import *

import pdb


class System():
    """System container class
    """
    def __init__(self, box=Box(None)):
        """
        """
        self.name = ''

        # numbas
        self.resids = list()
        self.resnames = np.empty(shape=0, dtype='str')
        self.n_molecules = 0

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

    def append_gbbs(self, gbbs, update_FF=True):
        """
        """
        for n_mols, gbb in enumerate(gbbs):
            n_new_atoms = gbb.types.shape[0]
            assert n_new_atoms == gbb.xyz.shape[0]
            if gbb.masses.shape[0] > 0:
                assert n_new_atoms == gbb.masses.shape[0]
            if gbb.charges.shape[0] > 0:
                assert n_new_atoms == gbb.charges.shape[0]

            # gather all the relevant numbers
            n_sys_atoms = self.xyz.shape[0]
            if self.types.shape[0] > 0:
                n_atom_types = max(self.pair_types.keys())
            else:
                n_atom_types = 0

            if self.bonds.shape[0] > 0:
                #n_bond_types = np.max(self.bonds[:, 0])
                n_bond_types = max(self.bond_types.keys())
            else:
                n_bond_types = 0

            if self.angles.shape[0] > 0:
                #n_angle_types = np.max(self.angles[:, 0])
                n_angle_types = max(self.angle_types.keys())
            else:
                n_angle_types = 0

            if self.dihedrals.shape[0] > 0:
                #n_dihedral_types = np.max(self.dihedrals[:, 0])
                n_dihedral_types = max(self.dihedral_types.keys())
            else:
                n_dihedral_types = 0

            if self.impropers.shape[0] > 0:
                #n_improper_types = np.max(self.impropers[:, 0])
                n_improper_types = max(self.improper_types.keys())
            else:
                n_improper_types = 0

            # update FF info for first molecule
            if (n_mols == 0) and (update_FF == True):
                # ---update numbering of pair coeffs
                convert_atoms = dict()
                for key, value in gbb.pair_types.items():
                    new_key = key + n_atom_types
                    convert_atoms[key] = new_key
                    self.pair_types[new_key] = value
                # ---update numbering of bond coeffs...
                convert_bonds = dict()
                for key, value in gbb.bond_types.items():
                    new_key = key + n_bond_types
                    convert_bonds[key] = new_key
                    self.bond_types[new_key] = value
                # ---update numbering of angle coeffs...
                convert_angles = dict()
                for key, value in gbb.angle_types.items():
                    new_key = key + n_angle_types
                    convert_angles[key] = new_key
                    self.angle_types[new_key] = value
                # ---update numbering of dihedral coeffs...
                convert_dihedrals = dict()
                for key, value in gbb.dihedral_types.items():
                    new_key = key + n_dihedral_types
                    convert_dihedrals[key] = new_key
                    self.dihedral_types[new_key] = value
                # ---update numbering of improper coeffs...
                convert_impropers = dict()
                for key, value in gbb.improper_types.items():
                    new_key = key + n_improper_types
                    convert_impropers[key] = new_key
                    self.improper_types[new_key] = value

            # ...and atom types
            for i, atom_type in enumerate(gbb.types):
                gbb.types[i] = convert_atoms[atom_type]

            # ...and bonds
            gbb.bonds[:, 1:] += n_sys_atoms
            for i, bond_type in enumerate(gbb.bonds[:, 0]):
                try:
                    gbb.bonds[i, 0] = convert_bonds[bond_type]
                except:
                    pdb.set_trace()

            # ...and angles
            gbb.angles[:, 1:] += n_sys_atoms
            for i, angle_type in enumerate(gbb.angles[:, 0]):
                gbb.angles[i, 0] = convert_angles[angle_type]

            # ...and dihedrals
            gbb.dihedrals[:, 1:] += n_sys_atoms
            for i, dihedral_type in enumerate(gbb.dihedrals[:, 0]):
                gbb.dihedrals[i, 0] = convert_dihedrals[dihedral_type]

            # ...and impropers
            gbb.impropers[:, 1:] += n_sys_atoms
            for i, improper_type in enumerate(gbb.impropers[:, 0]):
                gbb.impropers[i, 0] = convert_impropers[improper_type]

            self.types = np.append(self.types, gbb.types)
            self.bonds = np.vstack((self.bonds, gbb.bonds))
            self.angles = np.vstack((self.angles, gbb.angles))
            self.dihedrals = np.vstack((self.dihedrals, gbb.dihedrals))
            self.impropers = np.vstack((self.impropers, gbb.impropers))

            self.resids += [self.n_molecules + 1] * gbb.xyz.shape[0]

            self.masses = np.append(self.masses, gbb.masses)
            self.charges = np.append(self.charges, gbb.charges)
            self.xyz = np.vstack((self.xyz, gbb.xyz))
            self.n_molecules += 1

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
