import numpy as np
import re
import pdb


def read_frame_lammpstrj(trj, read_velocities=False):
    """Load a frame from a lammps dump file.

    Args:
        trj (file): LAMMPS dump file of format 'ID type x y z' or
                                               'ID type x y z vx vy vz'
        read_velocities (bool): if True, reads velocity data from file

    Returns:
        xyz (numpy.ndarray):
        types (numpy.ndarray):
        step (int):
        box (numpy.ndarray):
        vxyz (numpy.ndarray):
    """
    box = np.empty(shape=(3, 2))

    # --- begin header ---
    trj.readline()  # text
    step = int(trj.readline())  # timestep
    trj.readline()  # text
    n_atoms = int(trj.readline())  # num atoms
    trj.readline()  # text
    box[0] = trj.readline().split()  # x-dim of box
    box[1] = trj.readline().split()  # y-dim of box
    box[2] = trj.readline().split()  # z-dim of box
    trj.readline()  # text
    # --- end header ---

    xyz = np.empty(shape=(n_atoms, 3))
    xyz[:] = np.NAN
    types = np.empty(shape=(n_atoms), dtype='int')
    if read_velocities:
        vxyz = np.empty(shape=(n_atoms, 3))
        vxyz[:] = np.NAN

    # --- begin body ---
    for i in range(n_atoms):
        temp = trj.readline().split()
        a_ID = int(temp[0])  # atom ID
        types[a_ID - 1] = int(temp[1])  # atom type
        xyz[a_ID - 1] = map(float, temp[2:5])  # coordinates
        if read_velocities:
            vxyz[a_ID - 1] = map(float, temp[5:8])  # velocities
    # --- end body ---

    if read_velocities:
        return xyz, types, step, box, vxyz
    else:
        return xyz, types, step, box


def read_xyz(file_name):
    """Load an xyz file into a coordinate and a type array.
    """
    with open(file_name, 'r') as f:
        n_atoms = int(f.readline())  # num atoms
        f.readline()  # discard comment line

        xyz = np.empty(shape=(n_atoms, 3))
        types = np.empty(shape=(n_atoms), dtype='object')
        for i in range(n_atoms):
            temp = f.readline().split()
            types[i] = temp[0]  # atom type
            xyz[i] = map(float, temp[1:4])  # coordinates
        return xyz, types


def write_xyz(file_name, xyz, types, comment=''):
    """Write an xyz file.
    """
    assert xyz.shape[0] == types.shape[0]

    with open(file_name, 'w') as f:
        f.write(str(xyz.shape[0]) + '\n')  # num atoms
        f.write(comment + '\n')  # comment line

        for i, atom in enumerate(types):
            #        type   x     y     z
            f.write('%s %8.3f %8.3f %8.3f\n' %
                    (atom, xyz[i, 0], xyz[i, 1], xyz[i, 2]))
    print "Wrote file '" + file_name + "'"


def read_lammps_data(data_file, verbose=False):
    """Reads a LAMMPS data file

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
        -handling for comments
        -handling for directives not delimited by blank lines
        -allow specification of forcefield styles

    Args:
        data_file (str): name of lammps data file to read in
    Returns:
        lmp_data (dict):
            'xyz': xyz (numpy.ndarray)
            'types': types (numpy.ndarray)
            'masses': masses (numpy.ndarray)
            'charges': charges (numpy.ndarray)
            'bonds': bonds (numpy.ndarray)
            'angles': angles (numpy.ndarray)
            'dihedrals': dihedrals (numpy.ndarray)
            'pair_types': pair_types (dict)
            'bond_types': bond_types (dict)
            'angle_types': angle_types (dict)
            'dihedral_types': dihedral_type (dict)

        box (numpy.ndarray): box dimensions
    """
    bonds = np.empty(shape=(0, 3))
    angles = np.empty(shape=(0, 4))
    dihedrals = np.empty(shape=(0, 5))

    pair_types = dict()
    bond_types = dict()
    angle_types = dict()
    dihedral_types = dict()

    print "Reading '" + data_file + "'"
    with open(data_file, 'r') as f:
        data_lines = f.readlines()

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
                n_atoms = int(fields[0])
                xyz = np.empty(shape=(n_atoms, 3))
                types = np.empty(shape=(n_atoms), dtype='int')
                masses = np.empty(shape=(n_atoms))
                charges = np.empty(shape=(n_atoms))

            elif match.group('n_bonds'):
                fields = data_lines.pop(i).split()
                bonds = np.empty(shape=(float(fields[0]), 3), dtype='int')

            elif match.group('n_angles'):
                fields = data_lines.pop(i).split()
                angles = np.empty(shape=(float(fields[0]), 4), dtype='int')

            elif match.group('n_dihedrals'):
                fields = data_lines.pop(i).split()
                dihedrals = np.empty(shape=(float(fields[0]), 5), dtype='int')

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

                mass_dict = dict()  # type:mass

                #     not end of file         not blank line
                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    mass_dict[int(fields[0])] = float(fields[1])

            elif match.group('Atoms'):
                if verbose:
                    print 'Parsing Atoms...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    a_id = int(fields[0])
                    types[a_id - 1] = int(fields[2])
                    masses[a_id - 1] = mass_dict[int(fields[2])]
                    charges[a_id - 1] = float(fields[3])
                    xyz[a_id - 1] = np.array([float(fields[4]),
                                         float(fields[5]),
                                         float(fields[6])])

            elif match.group('PairCoeffs'):
                if verbose:
                    print 'Parsing Pair Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = data_lines.pop(i).split()
                    pair_types[int(fields[0])] = (float(fields[1]),
                                                  float(fields[2]))
            elif match.group('BondCoeffs'):
                if verbose:
                    print 'Parsing Bond Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(float, data_lines.pop(i).split())
                    bond_types[int(fields[0])] = fields[1:]

            elif match.group('AngleCoeffs'):
                if verbose:
                    print 'Parsing Angle Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(float, data_lines.pop(i).split())
                    angle_types[int(fields[0])] = fields[1:]

            elif match.group('DihedralCoeffs'):
                if verbose:
                    print 'Parsing Dihedral Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(float, data_lines.pop(i).split())
                    dihedral_types[int(fields[0])] = fields[1:]

            elif match.group('Bonds'):
                if verbose:
                    print 'Parsing Bonds...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(int, data_lines.pop(i).split())
                    bonds[fields[0] - 1] = fields[1:]

            elif match.group('Angles'):
                if verbose:
                    print 'Parsing Angles...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(int, data_lines.pop(i).split())
                    angles[fields[0] - 1] = fields[1:]

            elif match.group('Dihedrals'):
                if verbose:
                    print 'Parsing Dihedrals...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(int, data_lines.pop(i).split())
                    dihedrals[fields[0] - 1] = fields[1:]

            else:
                i += 1
        else:
            i += 1

    lmp_data = {'xyz': xyz,
                'types': types,
                'masses': masses,
                'charges': charges,
                'bonds': bonds,
                'angles': angles,
                'dihedrals': dihedrals,
                'pair_types': pair_types,
                'bond_types': bond_types,
                'angle_types': angle_types,
                'dihedral_types': dihedral_types
                }
    return lmp_data, box


def write_lammps_data(gbb, box, file_name='data.system', sys_name='system'):
    """Write gbb to lammps data file

    Args:
        gbb (Gbb): gbb object to write
        box (numpy.ndarray): box dimensions
        file_name (str): name of output file
        sys_name (str): name printed at top of data file
    """
    with open(file_name, 'w') as f:
        f.write(sys_name + '\n')
        f.write('\n')

        n_bonds = gbb.bonds.shape[0]
        n_angles = gbb.angles.shape[0]
        n_dihedrals = gbb.dihedrals.shape[0]

        f.write(str(gbb.xyz.shape[0]) + ' atoms\n')
        f.write(str(n_bonds) + ' bonds\n')
        f.write(str(n_angles) + ' angles\n')
        f.write(str(n_dihedrals) + ' dihedrals\n')
        f.write('\n')

        f.write(str(gbb.types.max()) + ' atom types\n')
        if n_bonds > 0:
            f.write(str(len(gbb.bond_types)) + ' bond types\n')
        if n_angles > 0:
            f.write(str(len(gbb.angle_types)) + ' angle types\n')
        if n_dihedrals > 0:
            f.write(str(len(gbb.dihedral_types)) + ' dihedral types\n')
        f.write('\n')

        f.write(str(box[0, 0]) + ' ' + str(box[0, 1]) + ' xlo xhi\n')
        f.write(str(box[1, 0]) + ' ' + str(box[1, 1]) + ' ylo yhi\n')
        f.write(str(box[2, 0]) + ' ' + str(box[2, 1]) + ' zlo zhi\n')

        f.write('\n')
        f.write('Masses\n')
        f.write('\n')

        # find unique masses and corresponding atomtypes
        masses = set()
        for i, mass in enumerate(gbb.masses):
            masses.add((gbb.types[i], mass))
        for mass in sorted(masses):
            f.write(" ".join(map(str, mass)) + '\n')

        f.write('\n')
        f.write('Atoms\n')
        f.write('\n')
        for i, coord in enumerate(gbb.xyz):
            f.write('%-6d %-6d %-6d %5.3f %8.3f %8.3f %8.3f\n'
                % (i+1,
                   gbb.mol_id,
                   gbb.types[i],
                   gbb.charges[i],
                   coord[0],
                   coord[1],
                   coord[2]))

        if n_bonds > 0:
            f.write('\n')
            f.write('Bonds\n')
            f.write('\n')
            for i, bond in enumerate(gbb.bonds):
                f.write(str(i+1) + " " + " ".join(map(str, bond)) + '\n')
        
        if n_angles > 0:
            f.write('\n')
            f.write('Angles\n')
            f.write('\n')
            for i, angle in enumerate(gbb.angles):
                f.write(str(i+1) + " " + " ".join(map(str, angle)) + '\n')

        if n_dihedrals > 0:
            f.write('\n')
            f.write('Dihedrals\n')
            f.write('\n')
            for i, dihedral in enumerate(gbb.dihedrals):
                f.write(str(i+1) + " " + " ".join(map(str, dihedral)) + '\n')

    print "Wrote file '" + file_name + "'"

