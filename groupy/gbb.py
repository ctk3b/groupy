import numpy as np
import copy
import re

from groupy.mdio import *
from groupy.general import *


class Gbb():
    """Generic building block.
    """
    def __init__(self):
        """
        """
        # numbas
        self.n_atoms = None
        self.n_bonds = None
        self.n_angles = None
        self.n_dihedrals = None

        self.mol_id = None

        # per atom
        self.types = None
        self.masses = None
        self.charges = None
        self.xyz = None

        # connectivity
        self.bonds = None
        self.angles = None
        self.dihedrals = None

        # forcefield info
        self.pair_types = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.dihedral_types = dict()

        # properties
        self.com = None

    # --- calculable properties ---
    def calc_com(self):
        """Calculate center of mass of a set of points.
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
        """Unwrap periodic boundary conditons of an object.

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

    def load_xyz(self, file_name):
        self.xyz, self.types = read_xyz(file_name)
        self.n_atoms = self.xyz.shape[0]

    def load_data(self, data_file, name='', verbose=False):
        """Reads a LAMMPS data file into a gbb object.

        *** Only works for directives delimited by blank lines ***

        Currently supports the following directives:
            Masses
            Pair Coeffs (must be mix geometry)
            Bond Coeffs (must be harmonic)
            Angle Coeffs (must be harmonic)
            Dihedral Coeffs (must be OPLS)
            Atoms
            Bonds
            Angles
            Dihedrals

        TODO:
            -add this function to mdio.py and add a call to populate gbb
            -handling for comments
            -handling for directives not delimited by blank lines
            -allow specification of forcefield styles

        Args:
            data_file (str): name of lammps data file to read in
            name (str): name for system
        Returns:
            box (numpy.ndarray): box dimensions
        """
        print "Reading '" + data_file + "'"
        with open(data_file, 'r') as f:
            data_lines = f.readlines()

        if name:
            self.name = name
        else:
            # toss out 'data.'
            prefix = re.search(r'\..+', data_file).group().strip('.')
            self.name = prefix

        directives = re.compile(r"""
            ((?P<n_atoms>.*atoms)
            |
            (?P<n_bonds>.*bonds)
            |
            (?P<n_angles>.*angles)
            |
            (?P<n_dihedrals>.*dihedrals)
            |
            (?P<box>.*xlo)
            |
            (?P<Masses>Masses)
            |
            (?P<PairCoeffs>Pair\sCoeffs)
            |
            (?P<BondCoeffs>Bond\sCoeffs)
            |
            (?P<AngleCoeffs>Angle\sCoeffs)
            |
            (?P<DihedralCoeffs>Dihedral\sCoeffs)
            |
            (?P<Atoms>Atoms)
            |
            (?P<Bonds>Bonds)
            |
            (?P<Angles>Angles)
            |
            (?P<Dihedrals>Dihedrals))
            """, re.VERBOSE)

        i = 0
        while i < len(data_lines):
            match = directives.match(data_lines[i])
            if match:
                if verbose:
                    print(match.groups())

                elif match.group('n_atoms'):
                    fields = data_lines.pop(i).split()
                    self.n_atoms = int(fields[0])
                    self.xyz = np.empty(shape=(self.n_atoms, 3))
                    self.types = np.empty(shape=(self.n_atoms))
                    self.masses = np.empty(shape=(self.n_atoms))
                    self.charges = np.empty(shape=(self.n_atoms))

                elif match.group('n_bonds'):
                    fields = data_lines.pop(i).split()
                    self.bonds = np.empty(shape=(float(fields[0]), 3))

                elif match.group('n_angles'):
                    fields = data_lines.pop(i).split()
                    self.angles = np.empty(shape=(float(fields[0]), 4))

                elif match.group('n_dihedrals'):
                    fields = data_lines.pop(i).split()
                    self.dihedrals = np.empty(shape=(float(fields[0]), 5))

                elif match.group('box'):
                    box = np.zeros(shape=(3, 2))
                    for j in range(3):
                        fields = map(float, data_lines.pop(i).split()[:2])
                        box[j, 0] = fields[0]
                        box[j, 1] = fields[1]

                elif match.group('Masses'):
                    if verbose:
                        print 'Parsing Masses...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    masses = dict()  # type:mass

                    #     not end of file         not blank line
                    while i < len(data_lines) and data_lines[i].strip():
                        fields = data_lines.pop(i).split()
                        masses[int(fields[0])] = float(fields[1])

                elif match.group('Atoms'):
                    if verbose:
                        print 'Parsing Atoms...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = data_lines.pop(i).split()
                        a_id = int(fields[0])
                        self.types[a_id - 1] = int(fields[2])
                        self.masses[a_id - 1] = masses[int(fields[2])]
                        self.charges[a_id - 1] = float(fields[3])
                        self.xyz[a_id - 1] = np.array([float(fields[4]),
                                             float(fields[5]),
                                             float(fields[6])])

                elif match.group('PairCoeffs'):
                    if verbose:
                        print 'Parsing Pair Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = data_lines.pop(i).split()
                        self.pair_types[int(fields[0])] = (float(fields[1]),
                                                      float(fields[2]))
                elif match.group('BondCoeffs'):
                    if verbose:
                        print 'Parsing Bond Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(float, data_lines.pop(i).split())
                        self.bond_types[int(fields[0])] = fields[1:]

                elif match.group('AngleCoeffs'):
                    if verbose:
                        print 'Parsing Angle Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(float, data_lines.pop(i).split())
                        self.angle_types[int(fields[0])] = fields[1:]

                elif match.group('DihedralCoeffs'):
                    if verbose:
                        print 'Parsing Dihedral Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(float, data_lines.pop(i).split())
                        self.dihedral_types[int(fields[0])] = fields[1:]

                elif match.group('Bonds'):
                    if verbose:
                        print 'Parsing Bonds...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(int, data_lines.pop(i).split())
                        self.bonds[fields[0] - 1] = fields[1:]

                elif match.group('Angles'):
                    if verbose:
                        print 'Parsing Angles...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(int, data_lines.pop(i).split())
                        self.angles[fields[0] - 1] = fields[1:]

                elif match.group('Dihedrals'):
                    if verbose:
                        print 'Parsing Dihedrals...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(int, data_lines.pop(i).split())
                        self.dihedrals[fields[0] - 1] = fields[1:]

                else:
                    i += 1
            else:
                i += 1
        return box
