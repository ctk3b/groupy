import numpy as np

from groupy.mdio import *
from groupy.general import *


class Gbb():
    """Generic building block.
    """
    def __init__(self):
        """
        """
        self.name = ''

        # numbas
        self.n_atoms = int()
        self.n_bonds = int()
        self.n_angles = int()
        self.n_dihedrals = int()
        self.mol_id = int()

        # per atom
        self.types = np.empty(shape=(0, 1), dtype='int')
        self.masses = np.empty(shape=(0, 1))
        self.charges = np.empty(shape=(0, 1))
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
        self.improper_types = dict()

        # properties
        self.com = float()
        self.r_gyr_sq = float()

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
    def translate(self, x=0.0, y=0.0, z=0.0, slow=False):
        """Translate by scalars
        """
        self.xyz[:, 0] += x
        self.xyz[:, 1] += y
        self.xyz[:, 2] += z

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
                        box_length = np.diff(box[k])
                        self.xyz[i, k] -= box_length * anint(dr / box_length)

    def wrap(self, box, dim=[True, True, True]):
        """Wrap coordinates for PBC.
        """

        for i, coords in enumerate(self.xyz):
            for k, c in enumerate(coords):  # TODO: vectorize
                if dim[k] == True:
                    if c < box[k, 0]:
                        self.xyz[i, k] = box[k, 1] - abs(box[k, 0] - c)
                    elif c > box[k, 1]:
                        self.xyz[i, k] = box[k, 0] + abs(c - box[k, 1])

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
    def load_mass(self, file_name):
        self.masses = np.loadtxt(file_name)

    def load_bonds(self, file_name):
        self.bonds = np.loadtxt(file_name)

    def load_xyz(self, file_name):
        self.xyz, self.types = read_xyz(file_name)
        self.n_atoms = self.xyz.shape[0]

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
