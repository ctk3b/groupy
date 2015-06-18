import collections
from copy import deepcopy

import numpy as np
from scipy.spatial import cKDTree

from box import *
from gbb import Gbb
from mdio import *


class System():
    """System container class - used for printing data files
    """
    def __init__(self, box=Box(None), gbbs=None, system_info=None):
        """
        """
        self.name = ''

        # individual components
        self.gbbs = list()

        # numbas
        self.resids = list()   # what is this
        self.resnames = np.empty(shape=0, dtype='str')

        # per atom properties
        """
        self.types = np.empty(shape=(0), dtype='int')
        self.masses = np.empty(shape=(0))
        self.charges = np.empty(shape=(0))
        self.xyz = np.empty(shape=(0, 3))
        """
        self.types = list()
        self.masses = list()
        self.charges = list()
        self.xyz = list()

        # connectivity
        """
        self.bonds = np.empty(shape=(0, 3), dtype='int')
        self.angles = np.empty(shape=(0, 4), dtype='int')
        self.dihedrals = np.empty(shape=(0, 5), dtype='int')
        self.impropers = np.empty(shape=(0, 5), dtype='int')
        """
        self.bonds = list()
        self.angles = list()
        self.dihedrals = list()
        self.impropers = list()

        # forcefield info 
        self.pair_types = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.dihedral_types = dict()

        # box
        self.box = box

        # self.system_info = list to give information on make up of system
        self._n_components = 0
        self.system_info = list()
        if system_info:
            for component in system_info:
                self.system_info.append(component)
                self._n_components += 1

        if gbbs:
            self.append_gbbs(gbbs)

        # other properties
        self.atom_offset = 0
        self.type_mass = dict()

    def n_atoms(self):
        n_atom = 0
        for component in self.system_info:
            n_atom += component[0] * component[1]
        return n_atom

    def n_molecules(self):
        n_molecule = 0
        for component in self.system_info:
            n_molecule += component[0]
        return n_molecule

    def add_components_to_info(self, t_components):
        for t_component in t_components:
            assert len(t_component) == 3
            self.system_info.append(t_component)
            self._n_components += 1

    def cumulative_atoms(self, n):
        """Return the number of atoms in the first n components.

        Args:
            n (int): return number of atoms in first n components
        """
        cum_atoms = 0
        for i in range(n):
            cum_atoms += self.system_info[i][0] * self.system_info[i][1]
        return cum_atoms

    def cumulative_molecules(self, n):
        """Return the number of molecules in the first n components.

        Args:
            n (int): return number of molecules in first n components
        """
        cum_molecules = 0
        for i in range(n):
            cum_molecules += self.system_info[i][0]
        return cum_molecules

    def print_system_info(self, filename="system-composition.txt"):
        f = open(filename, 'w')
        f.write('SYSTEM COMPOSITION\n\n')
        f.write('Species     Molecules     Atoms per molecule\n')
        f.write('-------     ---------     ------------------\n')
        for i, component in enumerate(self.system_info):
            f.write('%-12s%-14d%-d\n' % (component[2], component[0], component[1]))
        f.write('\n')
        f.write('%d total molecules\n' % self.n_molecules())
        f.write('%d total atoms' % self.n_atoms())
    def append_gbbs(self, gbbs, update_FF=True):
        """
            Append gbbs to the system, but don't enumerate yet

            This is a very simple function, unnecessary even?
            Args:
                gbbs: list of gbbs to add to system
        """
        for n_mols, gbb in enumerate(gbbs):
            # make sure there are the right number of atom properties
            n_new_atoms = len(gbb.types)
            assert n_new_atoms == len(gbb.xyz)
            if len(gbb.masses) > 0:
                assert n_new_atoms == len(gbb.masses)
            if len(gbb.charges) > 0:
                assert n_new_atoms == len(gbb.charges)

            self.gbbs.append(gbb)
            """
            if self.types.shape[0] > 0:
                n_atom_types = max(self.pair_types.keys())
            else:
                n_atom_types = 0
            """

            """
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
            """

            """
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

            # what is this?
            self.resids += [self.n_molecules + 1] * gbb.xyz.shape[0]

            self.masses = np.append(self.masses, gbb.masses)
            self.charges = np.append(self.charges, gbb.charges)
            self.xyz = np.vstack((self.xyz, gbb.xyz))
            self.n_molecules += 1
            """

    def enumerate_topology(self, destructive=True, int_types=False):
        # start from 0?
        if destructive == True:
            self.atom_offset = 0

        # add offset to each atom that gets printed
        import pdb
        for gbb in self.gbbs:
            for bond in gbb.bonds:
                for i, atom in enumerate(bond[:-1]):
                    bond[i+1] += self.atom_offset
                self.bonds.append(bond)
            for angle in gbb.angles:
                for i, atom in enumerate(angle[:-1]):
                    angle[i+1] += self.atom_offset
                self.angles.append(angle)
            for dihedral in gbb.dihedrals:
                for i, atom in enumerate(dihedral[:-1]):
                    dihedral[i+1] += self.atom_offset
                self.dihedrals.append(dihedral)
            for improper in gbb.impropers:
                for i, atom in enumerate(improper[:-1]):
                    improper[i+1] += self.atom_offset
                self.impropers.append(improper)
            self.atom_offset += gbb.xyz.shape[0]

            # oh and duh, add the atom-specific properties
            for pos in gbb.xyz:
                self.xyz.append(pos)
            for mass in gbb.masses:
                self.masses.append(mass)
            for atype in gbb.types:
                if not int_types:
                    self.types.append(atype)
                else:
                    self.types.append(int(atype))
            for charge in gbb.charges:
                self.charges.append(charge)

    def find_number_of_types(self):
        """Find unique atom, bond, etc... types.

        This could probably go along with enumerate_topology()
        Also creates a dict of {type: mass}
        """
        self.unique_atom_types = list()
        self.unique_bond_types = list()
        self.unique_angle_types = list()
        self.unique_dihedral_types = list()
        self.unique_improper_types = list()
        for i, atype in enumerate(self.types):
            if not atype in self.unique_atom_types:
                self.unique_atom_types.append(atype)
                self.type_mass[atype] = self.masses[i]
        for i, bond in enumerate(self.bonds):
            if not bond[0] in self.unique_bond_types:
                self.unique_bond_types.append(bond[0])
        for i, angle in enumerate(self.angles):
            if not angle[0] in self.unique_angle_types:
                self.unique_angle_types.append(angle[0])
        for i, dihedral in enumerate(self.dihedrals):
            if not dihedral[0] in self.unique_dihedral_types:
                self.unique_dihedral_types.append(dihedral[0])
        for i, improper in enumerate(self.impropers):
            if not improper[0] in self.unique_improper_types:
                self.unique_improper_types.append(improper[0])

        

    def print_lammpsdata(self, destructive=True, atom_style='full', 
            sys_name=None, filename='system.lammpsdata', ff_param_set=None,
            system_info=None):
        from groupy.mdio import write_lammpsdata
        self.resids = []
        for i, gbb in enumerate(self.gbbs):
            for atom in gbb.xyz:
                self.resids.append(i+1)
        self.enumerate_topology(destructive=destructive, int_types=True)
        self.find_number_of_types()
        write_lammpsdata(self, box=self.box, atom_style=atom_style,
                sys_name=sys_name, filename=filename, ff_param_set=ff_param_set,
                system_info=system_info)


    def write_lammpstrj(self, step=1, fmt='5col', 
            filename='traj.lammpstrj', mode='a'):
        self.enumerate_topology()
        write_lammpstrj_frame(np.asarray(self.xyz), np.asarray(self.types),
                step=1, box=self.box, fmt=fmt, filename=filename, mode=mode)


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
    
    def sort_by_name(self, order):
        """Sort the gbbs by name.

        Args:
            order: list of names by which to sort the gbbs.

        Useful for sorting molecules into blocks of gbbs with the same name.
        For example, you may want a mixed lipid system to be sorted so that lipi
        ds of type1 are the first N atoms, then type2 are the next M atoms and 
        so on and so forth.
        """
        self.sorted_gbbs = list()
        for lipid_name in order:
            for molecule in self.gbbs:
                if molecule.name == lipid_name:
                    self.sorted_gbbs.append(molecule)
        if len(self.sorted_gbbs) == len(self.gbbs):
            self.gbbs = self.sorted_gbbs
        else: 
            print "Warning: Not all molecules in sorted list, aborting sorting."

    def convert_from_traj(self, xyz, types, clear=True):
        """Convert a list of coordinates into a list of gbbs

        Args:
            xyz: numpy array of shape (N, 3) containing atomic positions
            types: numpy array of shape (N) containing atomic types
            clear: (if True) clear the system's gbb list
        """

        if np.all(clear):
            self.gbbs = list()
            
        # make sure xyz contains the right number of atoms and types
        error_message = 'Expected %d atoms, ' % self.n_atoms()
        error_message += 'received %d atoms.' % xyz.shape[0]
        assert xyz.shape[0] == self.n_atoms(), error_message
        error_message = '%d atoms given != ' % xyz.shape[0]
        error_message += '%d types given.' % types.shape[0]
        assert xyz.shape[0] == types.shape[0], error_message
        
        for i in range(self._n_components):
            for j in range(self.cumulative_atoms(i),
                           self.cumulative_atoms(i+1),
                           self.system_info[i][1]):
                t_gbb = Gbb(name=self.system_info[i][2])
                for k in range(self.system_info[i][1]):
                    t_gbb.insert_atom(pos=xyz[j+k], atype=types[j+k])
                self.gbbs.append(t_gbb)

    def get_composition(self):
        self.system_info = list()
        t_system_info = collections.OrderedDict()
        for i, gbb in enumerate(self.gbbs):
            if gbb.name in t_system_info:
                assert len(gbb.xyz) == t_system_info[gbb.name][1]
                t_system_info[gbb.name][0] += 1
            elif gbb.name not in t_system_info:
                t_system_info[gbb.name] = [1, len(gbb.xyz)]
        for comp in t_system_info:
            n_mol = t_system_info[comp][0]
            n_apm = t_system_info[comp][1]
            self.system_info.append([n_mol, n_apm, comp])

    def write_hoomd_xml(self, filename='start.xml'):
        self.enumerate_topology()
        write_hoomd_xml(self, self.box, filename=filename)
