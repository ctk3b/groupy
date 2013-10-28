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
        box_size (numpy.ndarray):
        vxyz (numpy.ndarray):
    """
    box_size = np.empty(shape=(3, 2))

    # --- begin header ---
    trj.readline()  # text
    step = int(trj.readline())  # timestep
    trj.readline()  # text
    n_atoms = int(trj.readline())  # num atoms
    trj.readline()  # text
    box_size[0] = trj.readline().split()  # x-dim of box
    box_size[1] = trj.readline().split()  # y-dim of box
    box_size[2] = trj.readline().split()  # z-dim of box
    trj.readline()  # text
    # --- end header ---

    xyz = np.empty(shape=(n_atoms, 3))
    xyz[:] = np.NAN
    types = np.empty(shape=(n_atoms))
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
        return xyz, types, step, box_size, vxyz
    else:
        return xyz, types, step, box_size


def read_xyz(file_name):
    """Load an xyz file into an array
    """
    with open(file_name, 'r') as f:
        n_atoms = int(f.readline())  # num atoms
        f.readline()  # discard comment line

        xyz = np.empty(shape=(n_atoms, 3))
        types = np.empty(shape=(n_atoms), dtype='string')
        for i in range(n_atoms):
            temp = f.readline().split()
            types[i] = temp[0]  # atom type
            xyz[i] = map(float, temp[1:4])  # coordinates
        return xyz, types
