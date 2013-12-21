import numpy as np

from mdio import *
from general import *
import copy

class Gbb():
    """Generic building block.
    """
    def __init__(self):
        """
        """
        self.name = ''
        self.mol_id = 1

        # numbas
        self.n_atoms = int()
        self.n_bonds = int()
        self.n_angles = int()
        self.n_dihedrals = int()

        self.resids = np.empty(shape=0, dtype='u4')
        self.resnames = np.empty(shape=0, dtype='str')

        # per atom
        self.types = np.empty(shape=(0, 1))
        self.masses = np.empty(shape=(0, 1))
        self.charges = np.empty(shape=(0, 1))
        self.xyz = np.empty(shape=(0, 3))
        self.vel = np.empty(shape=(0, 3))

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
        self.improper_types = dict()

        # properties
        self.com = np.empty(shape=(3))
        self.r_gyr_sq = float()

    def delete_ions(self, ids):
        """Removes all information associated selected indices

        NOTE: does not work for particles bonded to anything
            
        """
        ids_1 = ids + 1  # 1-indexed list for atom numbering
        atom_map = dict()
        count = 1
        for i in range(1, self.xyz.shape[0] + 1):
            if i in ids_1:
                pass 
            else:
                atom_map[i] = count
                count += 1

        for bond in self.bonds:
            bond[1] = atom_map[bond[1]]
            bond[2] = atom_map[bond[2]]

        for angle in self.angles:
            angle[1] = atom_map[angle[1]]
            angle[2] = atom_map[angle[2]]
            angle[3] = atom_map[angle[3]]

        for dihedral in self.dihedrals:
            dihedral[1] = atom_map[dihedral[1]]
            dihedral[2] = atom_map[dihedral[2]]
            dihedral[3] = atom_map[dihedral[3]]
            dihedral[4] = atom_map[dihedral[4]]

        for improper in self.impropers:
            improper[1] = atom_map[improper[1]]
            improper[2] = atom_map[improper[2]]
            improper[3] = atom_map[improper[3]]
            improper[4] = atom_map[improper[4]]

        self.types = np.delete(self.types, ids)
        self.charges = np.delete(self.charges, ids)
        self.masses = np.delete(self.masses, ids)
        self.xyz = np.delete(self.xyz, ids, 0)
        if self.vel.shape[0] > 0:
            self.vel = np.delete(self.vel, ids, 0)


    def delete_residue(self, i):
        """Removes all information associated with residue 'i'

        NOTE: currently does not affect connectivity
        """
        keep_these = (self.resids != i)
        self.resids = self.resids[keep_these]
        self.resnames = self.resnames[keep_these]
        self.types = self.types[keep_these]
        self.xyz = self.xyz[keep_these]
        self.vel = self.vel[keep_these]

    # --- calculable properties ---
    def calc_com(self):
        """Calculate center of mass of a set of points.
        """
        assert self.xyz.shape[0] == self.masses.shape[0]
        cum_mass = 0
        com = np.array([0., 0., 0.])
        com = np.sum(self.xyz.T * self.masses, axis=1)
        cum_mass = np.sum(self.masses)
        com /= cum_mass
        self.com = com

    def calc_inertia_tensor(self):
        """Calculates the moment of inertia tensor of the gbb.

        Returns:
            I (np.ndarray): moment of inertia tensor
        """
        assert self.xyz.shape[0] == self.masses.shape[0]
        I = np.zeros(shape=(3, 3))
        self.calc_com()
        for i, coord0 in enumerate(self.xyz):
            mass = self.masses[i]
            coord = coord0 - self.com
            I[0, 0] += mass * (coord[1] * coord[1] + coord[2] * coord[2])
            I[1, 1] += mass * (coord[0] * coord[0] + coord[2] * coord[2])
            I[2, 2] += mass * (coord[0] * coord[0] + coord[1] * coord[1])
            I[0, 1] -= mass * coord[0] * coord[1]
            I[0, 2] -= mass * coord[0] * coord[2]
            I[1, 2] -= mass * coord[1] * coord[2]
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        return I

    def calc_r_gyr_sq(self):
        """Calculate radius of gyration.
        """
        pos_mean = np.mean(self.xyz, axis=0)
        r_gyr_sq = np.sum(np.sum(np.abs(self.xyz - pos_mean)**2, axis=1))
        r_gyr_sq /= self.xyz.shape[0]
        self.r_gyr_sq =  r_gyr_sq

    # --- manipulations ---
    def translate(self, xyz=np.zeros(3)):
        """Translate by scalars
        """
        self.xyz[:, 0] += xyz[0]
        self.xyz[:, 1] += xyz[1]
        self.xyz[:, 2] += xyz[2]

    def shift_com_to_origin(self):
        self.calc_com()
        self.translate(-self.com)

    def rotate(self, angles=[0.0, 0.0, 0.0]):
        """Rotate around given axes by given angles
        """
        if angles[0] != 0.0:
            theta = angles[0]
            T = np.eye(3)
            T[1, 1] = np.cos(theta)
            T[1, 2] = -np.sin(theta)
            T[2, 1] = np.sin(theta)
            T[2, 2] = np.cos(theta)
            self.xyz = np.dot(self.xyz, T.T)
        if angles[1] != 0.0:
            theta = angles[1]
            T = np.eye(3)
            T[0, 0] = np.cos(theta)
            T[0, 2] = np.sin(theta)
            T[2, 0] = -np.sin(theta)
            T[2, 2] = np.cos(theta)
            self.xyz = np.dot(self.xyz, T.T)
        if angles[2] != 0.0:
            theta = angles[2]
            T = np.eye(3)
            T[0, 0] = np.cos(theta)
            T[0, 1] = -np.sin(theta)
            T[1, 0] = np.sin(theta)
            T[1, 1] = np.cos(theta)
            self.xyz = np.dot(self.xyz, T.T)

    def scale():
        pass

    def unwrap(self, box, dim=[True, True, True]):
        """Unwrap periodic boundary conditons.

        Requires that object being unwrapped does not span more than half the
        box length.
        """
        for i, atom in enumerate(self.xyz):
            if i != 0:
                for k in range(3):  # TODO: vectorize
                    if dim[k] == True:
                        dr = self.xyz[i, k] - self.xyz[0, k]
                        box_length = box.length[k]
                        self.xyz[i, k] -= box_length * anint(dr / box_length)

    def wrap(self, box, dim=[True, True, True]):
        """Wrap coordinates for PBC.
        """

        for i, coords in enumerate(self.xyz):
            for k, c in enumerate(coords):  # TODO: vectorize
                if dim[k] == True:
                    if c < box.dims[k, 0]:
                        self.xyz[i, k] = box.dims[k, 1] - abs(box.dims[k, 0] - c)
                    elif c > box.dims[k, 1]:
                        self.xyz[i, k] = box.dims[k, 0] + abs(c - box.dims[k, 1])

    def mirror(self, box, dim=[False, False, True], d=0.0, verbose=False):
        """Creates a mirror image of gbb in a given direction.

        Intended to prepare SAM systems for shearing.

        *** currently only flips in the z-direction ***
        *** assume SAM is oriented upwards in the positive z direction ***
        *** currently bases separation distance on max z coordinate ***

        TODO:
            flip in any direction based on user input

        Args:
            d (float): desired separation distance
        """
        # find max/min z
        z_max = self.xyz[:, 2].max()
        z_min = self.xyz[:, 2].min()
        # adjust box size
        box[2] = [z_min - 1, z_max + abs(z_max - z_min) + d + 1]

        # duplicate types entries
        self.types = np.hstack((self.types, self.types))
        # make atom copies
        new_xyz = np.empty_like(self.xyz)
        for i, coord in enumerate(new_xyz):
            # adjust coords
            old = self.xyz[i]
            new_xyz[i, 0] = box[0, 0] + abs(box[0, 1] - old[0])
            new_xyz[i, 1] = box[1, 0] + abs(box[1, 1] - old[1])
            new_xyz[i, 2] = z_max + d + abs(z_max - old[2])
        self.xyz = np.vstack((self.xyz, new_xyz))
        """
        # make bond copies
        if self.bonds:
            n_bonds = len(self.bonds)
            for i, bond in enumerate(self.bonds):
                if i == n_bonds:
                    break
                new_bond = [len(self.bonds) + 1, bond[1], 0, 0]
                for j in range(2):
                    new_bond[j + 2] = bond[j + 2] + n_atoms
                self.add_bond(new_bond)

        # make angle copies
        if self.angles:
            n_angles = len(self.angles)
            for i, angle in enumerate(self.angles):
                if i == n_angles:
                    break
                new_angle = [len(self.angles) + 1, angle[1], 0, 0, 0]
                for j in range(3):
                    new_angle[j + 2] = angle[j + 2] + n_atoms
                self.add_angle(new_angle)

        # make dihedral copies
        if self.dihedrals:
            n_dihedrals = len(self.dihedrals)
            for i, dihedral in enumerate(self.dihedrals):
                if i == n_dihedrals:
                    break
                new_dihedral = [len(self.dihedrals) + 1, dihedral[1], 0, 0, 0, 0]
                for j in range(4):
                    new_dihedral[j + 2] = dihedral[j + 2] + n_atoms
                self.add_dihedral(new_dihedral)
        """

    # --- io ---
    def load_prototypes(self, base_name, prototypes):
        if 'mass' in prototypes:
            self.load_mass(base_name + '_mass.txt')
        if 'charge' in prototypes:
            self.load_charge(base_name + '_charge.txt')
        if 'coord' in prototypes:
            self.load_coord(base_name + '_coord.txt')
        if 'bond' in prototypes:
            self.load_bond(base_name + '_bond.txt')
        if 'angle' in prototypes:
            self.load_angle(base_name + '_angle.txt')
        if 'dihedral' in prototypes:
            self.load_dihedral(base_name + '_dihedral.txt')
        if 'pair_types' in prototypes:
            self.load_pair_types(base_name + '_pair_types.txt')
        if 'bond_types' in prototypes:
            self.load_bond_types(base_name + '_bond_types.txt')
        if 'angle_types' in prototypes:
            self.load_angle_types(base_name + '_angle_types.txt')
        if 'dihedral_types' in prototypes:
            self.load_dihedral_types(base_name + '_dihedral_types.txt')

    def load_mass(self, file_name):
        self.masses = np.loadtxt(file_name)

    def load_bond(self, file_name):
        self.bonds = np.loadtxt(file_name, dtype='int')

    def load_angle(self, file_name):
        self.angles = np.loadtxt(file_name, dtype='int')

    def load_dihedral(self, file_name):
        self.dihedrals = np.loadtxt(file_name, dtype='int')

    def load_improper(self, file_name):
        self.impropers = np.loadtxt(file_name, dtype='int')

    def load_pair_types(self, file_name):
        data = np.loadtxt(file_name)
        self.pair_types = dict(zip(data[:, 0].astype(int), 
            zip(data[:, 1], data[:, 2])))

    def load_bond_types(self, file_name):
        data = np.loadtxt(file_name)
        self.bond_types = dict(zip(data[:, 0].astype(int), 
            zip(data[:, 1], data[:, 2])))

    def load_angle_types(self, file_name):
        data = np.loadtxt(file_name)
        self.angle_types = dict(zip(data[:, 0].astype(int), 
            zip(data[:, 1], data[:, 2])))

    def load_dihedral_types(self, file_name):
        data = np.loadtxt(file_name)
        self.dihedral_types = dict(zip(data[:, 0].astype(int), 
            zip(data[:, 1], data[:, 2], data[:, 3], data[:, 4])))

    def load_charge(self, file_name):
        self.charges = np.loadtxt(file_name)

    def load_xyz(self, file_name):
        self.xyz, self.types = read_xyz(file_name)
        self.n_atoms = self.xyz.shape[0]

    def load_coord(self, file_name):
        coords = np.loadtxt(file_name)
        self.types = coords[:, 0].astype(int)
        self.xyz = coords[:, 1:]

    def load_gro(self, file_name):
        self.resids, self.resnames, self.types, self.xyz, self.vel, box= read_gro(file_name)
        self.n_atoms = self.xyz.shape[0]
        return box

    def load_lammps_data(self, data_file, verbose=False):
        lmp_data, box = read_lammps_data(data_file)
        self.xyz = lmp_data['xyz']
        self.types = lmp_data['types']
        self.masses = lmp_data['masses']
        self.charges = lmp_data['charges']
        self.bonds = lmp_data['bonds']
        self.angles = lmp_data['angles']
        self.dihedrals = lmp_data['dihedrals']
        self.pair_types = lmp_data['pair_types']
        self.bond_types = lmp_data['bond_types']
        self.angle_types = lmp_data['angle_types']
        self.dihedral_types = lmp_data['dihedral_types']
        return box
