import numpy as np

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
    n = int(trj.readline())  # num atoms
    trj.readline()  # text
    box_size[0] = trj.readline().split()  # x-dim of box
    box_size[1] = trj.readline().split()  # y-dim of box
    box_size[2] = trj.readline().split()  # z-dim of box
    trj.readline()  # text
    # --- end header ---

    xyz = np.empty(shape=(n, 3))
    xyz[:] = np.NAN
    types = np.empty(shape=(n))
    if read_velocities:
        vxyz = np.empty(shape=(n, 3))
        vxyz[:] = np.NAN

    # --- begin body ---
    for i in range(n):
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
