import numpy as np
import warnings
import re
import pdb

from box import *


def read_frame_lammpstrj(trj, read_velocities=False):
    """Load a frame from a LAMMPS dump file.

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
        return xyz, types, step, Box(box[:, 1], box[:, 0]), vxyz
    else:
        return xyz, types, step, Box(box[:, 1], box[:, 0])


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


def write_xyz(xyz, types, file_name, comment=''):
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
        data_file (str): name of LAMMPS data file to read in
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
    bonds = np.empty(shape=(0, 3), dtype='int')
    angles = np.empty(shape=(0, 4), dtype='int')
    dihedrals = np.empty(shape=(0, 5), dtype='int')
    impropers = np.empty(shape=(0, 5), dtype='int')

    pair_types = dict()
    bond_types = dict()
    angle_types = dict()
    dihedral_types = dict()
    improper_types = dict()

    print "Reading '" + data_file + "'"
    with open(data_file, 'r') as f:
        data_lines = f.readlines()

    # TODO: improve robustness of xlo regex
    directives = re.compile(r"""
        ((?P<n_atoms>\s*\d+\s+atoms)
        |
        (?P<n_bonds>\s*\d+\s+bonds)
        |
        (?P<n_angles>\s*\d+\s+angles)
        |
        (?P<n_dihedrals>\s*\d+\s+dihedrals)
        |
        (?P<box>.+xlo)
        |
        (?P<Masses>\s*Masses)
        |
        (?P<PairCoeffs>\s*Pair\sCoeffs)
        |
        (?P<BondCoeffs>\s*Bond\sCoeffs)
        |
        (?P<AngleCoeffs>\s*Angle\sCoeffs)
        |
        (?P<DihedralCoeffs>\s*Dihedral\sCoeffs)
        |
        (?P<ImproperCoeffs>\s*Improper\sCoeffs)
        |
        (?P<Atoms>\s*Atoms)
        |
        (?P<Bonds>\s*Bonds)
        |
        (?P<Angles>\s*Angles)
        |
        (?P<Dihedrals>\s*Dihedrals)
        |
        (?P<Impropers>\s*Impropers))
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
                dims = np.zeros(shape=(3, 2))
                for j in range(3):
                    fields = map(float, data_lines.pop(i).split()[:2])
                    dims[j, 0] = fields[0]
                    dims[j, 1] = fields[1]
                box = Box(dims[:, 1], dims[:, 0])

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
                    if len(fields) == 7:
                        a_id = int(fields[0])
                        types[a_id - 1] = int(fields[2])
                        masses[a_id - 1] = mass_dict[int(fields[2])]
                        charges[a_id - 1] = float(fields[3])
                        xyz[a_id - 1] = np.array([float(fields[4]),
                                             float(fields[5]),
                                             float(fields[6])])

                    # non-official file format
                    if len(fields) == 8:
                        a_id = int(fields[0])
                        types[a_id - 1] = int(fields[1])
                        masses[a_id - 1] = mass_dict[int(fields[1])]
                        charges[a_id - 1] = 0.0
                        xyz[a_id - 1] = np.array([float(fields[2]),
                                             float(fields[3]),
                                             float(fields[4])])


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

            elif match.group('ImproperCoeffs'):
                if verbose:
                    print 'Parsing Improper Coeffs...'
                data_lines.pop(i)
                data_lines.pop(i)

                while i < len(data_lines) and data_lines[i].strip():
                    fields = map(float, data_lines.pop(i).split())
                    improper_types[int(fields[0])] = fields[1:]

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

            elif match.group('Impropers'):
                if verbose:
                    print 'Parsing Impropers...'
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
                'impropers': impropers,
                'pair_types': pair_types,
                'bond_types': bond_types,
                'angle_types': angle_types,
                'dihedral_types': dihedral_types,
                'improper_types': improper_types
                }
    return lmp_data, box


def write_lammps_data(gbb, box=None, file_name='data.system', sys_name='system', prototype=False):
    """Write gbb to LAMMPS data file

    Args:
        gbb (Gbb): gbb object to write
        box (numpy.ndarray): box dimensions
        file_name (str): name of output file
        sys_name (str): name printed at top of data file
    """
    if (box == None) and np.sum(gbb.box.dims) > 0:
        box = gbb.box
    elif (box == None):
        raise Exception("Box is empty")

    with open(file_name, 'w') as f:
        f.write(sys_name + '\n')
        f.write('\n')

        n_bonds = int(gbb.bonds.shape[0])
        n_angles = int(gbb.angles.shape[0])
        n_dihedrals = int(gbb.dihedrals.shape[0])
        n_impropers = int(gbb.impropers.shape[0])

        f.write(str(gbb.xyz.shape[0]) + ' atoms\n')
        f.write(str(n_bonds) + ' bonds\n')
        f.write(str(n_angles) + ' angles\n')
        f.write(str(n_dihedrals) + ' dihedrals\n')
        f.write('\n')

        f.write(str(int(gbb.types.max())) + ' atom types\n')
        if n_bonds > 0:
            f.write(str(int(gbb.bonds[:, 0].max())) + ' bond types\n')
        if n_angles > 0:
            f.write(str(int(gbb.angles[:, 0].max())) + ' angle types\n')
        if n_dihedrals > 0:
            f.write(str(int(gbb.dihedrals[:, 0].max())) + ' dihedral types\n')
        f.write('\n')

        f.write('%-8.4f %8.4f xlo xhi\n' % (box.dims[0, 0], box.dims[0, 1]))
        f.write('%-8.4f %8.4f ylo yhi\n' % (box.dims[1, 0], box.dims[1, 1]))
        f.write('%-8.4f %8.4f zlo zhi\n' % (box.dims[2, 0], box.dims[2, 1]))

        f.write('\n')
        f.write('Masses\n')
        f.write('\n')

        # find unique masses and corresponding atomtypes
        masses = set()
        if prototype:
            f_mass = open(sys_name + '_mass.txt', 'w')
        for i, mass in enumerate(gbb.masses):
            if prototype:
                f_mass.write('%8.4f\n' % (mass))
            masses.add((int(gbb.types[i]), mass))
        for mass in sorted(masses):
            f.write(" ".join(map(str, mass)) + '\n')
        if prototype:
            f_mass.close()

        if len(gbb.pair_types) > 0:
            if prototype:
                f_pair_type = open(sys_name + '_pair_types.txt', 'w')
            f.write('\n')
            f.write('Pair Coeffs\n')
            f.write('\n')
            for i, value in sorted(gbb.pair_types.items()):
                f.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))
                if prototype:
                    f_pair_type.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))
            if prototype:
                f_pair_type.close()

        if len(gbb.bond_types) > 0:
            if prototype:
                f_bond_type = open(sys_name + '_bond_types.txt', 'w')
            f.write('\n')
            f.write('Bond Coeffs\n')
            f.write('\n')
            for i, value in sorted(gbb.bond_types.items()):
                f.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))
                if prototype:
                    f_bond_type.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))
            if prototype:
                f_bond_type.close()


        if len(gbb.angle_types) > 0:
            if prototype:
                f_ang_type = open(sys_name + '_angle_types.txt', 'w')
            f.write('\n')
            f.write('Angle Coeffs\n')
            f.write('\n')
            for i, value in sorted(gbb.angle_types.items()):
                f.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))
                if prototype:
                    f_ang_type.write("%d %8.4f %8.4f\n" % (i, value[0], value[1]))
            if prototype:
                f_ang_type.close()


        if len(gbb.dihedral_types) > 0:
            if prototype:
                f_dih_type = open(sys_name + '_dihedral_types.txt', 'w')
            f.write('\n')
            f.write('Dihedral Coeffs\n')
            f.write('\n')
            for i, value in sorted(gbb.dihedral_types.items()):
                f.write("%d %8.4f %8.4f %8.4f %8.4f\n" % (i, value[0], value[1], value[2], value[3]))
                if prototype:
                    f_dih_type.write("%d %8.4f %8.4f %8.4f %8.4f\n" % (i, value[0], value[1], value[2], value[3]))
            if prototype:
                f_dih_type.close()



        if prototype:
            f_coord = open(sys_name + '_coord.txt', 'w')
            f_charge = open(sys_name + '_charge.txt', 'w')
        f.write('\n')
        f.write('Atoms\n')
        f.write('\n')
        for i, coord in enumerate(gbb.xyz):
            if len(gbb.resids) > 0:
                resid = gbb.resids[i]
            elif len(gbb.resids) == 0:
                resid = 1
            f.write('%-6d %-6d %-6d %5.8f %8.3f %8.3f %8.3f\n'
                % (i+1,
                   resid,
                   gbb.types[i],
                   gbb.charges[i],
                   coord[0],
                   coord[1],
                   coord[2]))
            if prototype:
                f_coord.write('%d %8.4f %8.4f %8.4f\n' % (gbb.types[i],
                    coord[0], coord[1], coord[2]))
                f_charge.write('%8.4f\n' % (gbb.charges[i]))
        if prototype:
            f_coord.close()
            f_charge.close()

        if n_bonds > 0:
            if prototype:
                f_bond = open(sys_name + '_bond.txt', 'w')
            f.write('\n')
            f.write('Bonds\n')
            f.write('\n')
            for i, bond in enumerate(gbb.bonds):
                f.write(str(i+1) + " " + " ".join(map(str, bond)) + '\n')
                if prototype:
                    f_bond.write(' '.join(map(str, bond)) + '\n')
            if prototype:
                f_bond.close()

        if n_angles > 0:
            if prototype:
                f_angle = open(sys_name + '_angle.txt', 'w')
            f.write('\n')
            f.write('Angles\n')
            f.write('\n')
            for i, angle in enumerate(gbb.angles):
                f.write(str(i+1) + " " + " ".join(map(str, angle)) + '\n')
                if prototype:
                    f_angle.write(' '.join(map(str, angle)) + '\n')
            if prototype:
                f_angle.close()

        if n_dihedrals > 0:
            if prototype:
                f_dihedral = open(sys_name + '_dihedral.txt', 'w')
            f.write('\n')
            f.write('Dihedrals\n')
            f.write('\n')
            for i, dihedral in enumerate(gbb.dihedrals):
                f.write(str(i+1) + " " + " ".join(map(str, dihedral)) + '\n')
                if prototype:
                    f_dihedral.write(' '.join(map(str, dihedral)) + '\n')
            if prototype:
                f_dihedral.close()

        if n_impropers > 0:
            if prototype:
                f_dihedral = open(sys_name + '_dihedral.txt', 'w')
            f.write('\n')
            f.write('Impropers\n')
            f.write('\n')
            for i, improper in enumerate(gbb.impropers):
                f.write(str(i+1) + " " + " ".join(map(str, improper)) + '\n')
                if prototype:
                    f_dihedral.write(' '.join(map(str, dihedral)) + '\n')

            if prototype:
                f_dihedral.close()

    print "Wrote file '" + file_name + "'"

def read_gro(file_name):
    """
    """
    if not file_name.endswith('.gro'):
        warnings.warn("File name passed to read_gro() does not end with '.gro'")

    with open(file_name, 'r') as f:
        sys_name = f.readline().strip()
        n_atoms = int(f.readline())

        resids = np.empty(shape=(n_atoms), dtype='u4')
        resnames = np.empty(shape=(n_atoms), dtype='a5')
        types = np.empty(shape=(n_atoms), dtype='a5')
        xyz = np.empty(shape=(n_atoms, 3), dtype='f2')
        vel = np.empty(shape=(n_atoms, 3), dtype='f2')
        for i in range(n_atoms):
            line = f.readline()
            resids[i] = int(line[:5])
            resnames[i] = line[5:10].strip()
            types[i] = line[10:15].strip()
            xyz[i, 0] = float(line[20:28])
            xyz[i, 1] = float(line[28:36])
            xyz[i, 2] = float(line[36:44])
            try:
                vel[i, 0] = float(line[44:52])
                vel[i, 1] = float(line[52:60])
                vel[i, 2] = float(line[60:68])
            except:
                vel[i, 0] = 0.0
                vel[i, 1] = 0.0
                vel[i, 2] = 0.0

        box = np.zeros(shape=(3, 2))
        line = map(float, f.readline().split())
        if len(line) == 3:
            box[0, 1] = line[0]
            box[1, 1] = line[1]
            box[2, 1] = line[2]

    print "Read file '" + file_name + "'"
    return resids, resnames, types, xyz, vel, box




def write_gro(gbb, box, grofile='system.gro', sys_name='system'):
    """Write gbb to GROMACS .gro file

    Note: to my knowledge GROMACS only deals with positive coordinate
    values. Use method 'positive_coords' to shift coordinates appropriately

    Args:
        grofile (str): name of .gro structure file
    """
    with open(grofile, 'w') as f:
        f.write(sys_name + '\n')
        f.write(str(gbb.xyz.shape[0]) + '\n')
        for i, coord in enumerate(gbb.xyz):
            f.write('%5d%-4s%6s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n'
                      %(gbb.resids[i],
                        gbb.resnames[i],
                        gbb.types[i],
                        i+1,
                        coord[0],
                        coord[1],
                        coord[2],
                        gbb.vel[i, 0],
                        gbb.vel[i, 1],
                        gbb.vel[i, 2]))
        f.write('%10.5f%10.5f%10.5f\n'
               %(box[0, 1],
                 box[1, 1],
                 box[2, 1]))
    print "Wrote file '" + grofile + "'"


def write_top(gbb, topfile='system.top'):
    """Write forcefield information to GROMACS .top file

    *** THIS IS ONLY THE SKELETON FOR THE FORMAT ***

    TODO:
        -fill in body 
    """
    with open(topfile, 'w') as f:

        # defaults
        f.write('[ defaults ]\n')
        f.write('%d%16d%18s%20.4f%8.4f\n'
                %(1, 3, 'yes', 1, 1)) #  TODO: options need global storage 

        # atomtypes
        f.write('[ atomtypes ]\n')
        f.write('\n')

        # moleculetype
        for a in blah:  # BROKED
            f.write('[ moleculetype ]\n')
            f.write('\n')
            f.write('[ atoms ]\n')
            f.write('\n')

            f.write('[ bonds ]\n')
            for bond in gbb.bonds:
                r = gbb.bond_types[bond][1]
                k = gbb.bond_types[bond][2]
                f.write('%6d%7d%4d%18.8e%18.8e\n'
                        %(bond[2], bond[3], 1, r, k))

            f.write('[ angles ]\n')
            f.write('\n')
            f.write('[ dihedrals ]\n')
            f.write('\n')

        # system
        f.write('[ system ]\n')

        # molecules
        f.write('[ molecules ]\n')
    print "Wrote file '" + topfile + "'"
