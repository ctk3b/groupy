import numpy as np
import pdb


def read_frame_lammpstrj(trj, read_velocities=False):
    """Load a frame from a lammps dump file

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
    """Load an xyz file into a coordinate and a type array
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
    """Write an xyz file
    """
    assert xyz.shape[0] == types.shape[0]

    with open(file_name, 'w') as f:
        f.write(str(xyz.shape[0]) + '\n')  # num atoms
        f.write(comment + '\n')  # comment line

        for i, atom in enumerate(types):
            # type  x y z
            f.write('%s %8.3f %8.3f %8.3f\n' %  
                    (atom, xyz[i, 0], xyz[i, 1], xyz[i, 2]))
