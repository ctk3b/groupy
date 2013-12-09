import pdb

from groupy.mdio import *
from groupy.general import *


def calc_rdf(file_name, pairs=None, n_bins=100, max_frames=np.inf, opencl=False):
    """Radial distribution function - g(r)

    Args:
        file_name (str): name of trajectory file
        pairs (list): pair of types between which to calculate g(r)
        n_bins (int):
        max_frames (int):
    Returns:
        r (np.ndarray): radii values corresponding to bins
        g_r (np.ndarray): radial distribution functions at radii, r
    """
    if opencl:
        import pyopencl as cl

        ctx = cl.create_some_context()
        queue = cl.CommandQueue(ctx)
        mf = cl.mem_flags

        with open('g_r.cl', 'r') as f:
            source = "".join(f.readlines())
        program = cl.Program(ctx, source).build()

    r_range = (0.0, 8.0)
    g_r, edges = np.histogram([0], bins=n_bins, range=r_range)
    g_r[0] = 0
    g_r = g_r.astype(np.float64)
    total_volume = 0
    n_frames = 0
    rho = 0

    with open(file_name, 'r') as trj:
        while n_frames < max_frames:
            try:
                xyz, types, _, box = read_frame_lammpstrj(trj)
            except:
                print "Reach end of '" + file_name + "'"
                break
            n_frames += 1
            n_atoms = xyz.shape[0]
            print "read " + str(n_frames)

            if opencl:
                for pair in pairs:
                    x = np.asarray(xyz[:, 0], dtype='float32')
                    y = np.asarray(xyz[:, 1], dtype='float32')
                    z = np.asarray(xyz[:, 2], dtype='float32')
                    # in
                    x_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=x)
                    y_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=y)
                    z_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=z)
                    # out
                    g_buf = cl.Buffer(ctx, mf.WRITE_ONLY, x.nbytes * x.nbytes)

                    # run
                    program.g_r(queue,
                            (x.shape[0] * x.shape[0], ),
                            None,
                            x_buf,
                            y_buf,
                            z_buf,
                            l_buf,
                            d_buf)

                    # extract
                    g_buf = np.empty(x.shape[0] * x.shape[0], dtype=np.float32)
                    cl.enqueue_read_buffer(queue, g_buf, g).wait()

                    pdb.set_trace()

                    temp_g_r, _ = np.histogram(d, bins=n_bins, range=r_range)
                    g_r += temp_g_r

            # TODO: loop over multiple pairs for pure numpy version
            else:
                if pairs is None:
                #if pairs is 'asdfa':
                    d = np.zeros(n_atoms * n_atoms)
                    for i, xyz_i in enumerate(xyz):
                        xyz_j = np.vstack([xyz[:i], xyz[i+1:]])
                        temp_d = calc_distance_pbc(xyz_i, xyz_j, box.length)
                        start = (n_atoms - 1) * i
                        stop = (n_atoms - 1) * (i + 1)
                        d[start:stop] = temp_d
                    temp_g_r, _ = np.histogram(d, bins=n_bins, range=r_range)

                    g_r += temp_g_r
                    pdb.set_trace()
                    if n_frames == 1:
                        print g_r
                fast = temp_d
                # all-all
                #elif pairs is None:
                slow = temp
                if pairs is None:
                #elif pairs is 'asfd':
                    for i, xyz_i in enumerate(xyz):
                        xyz_j = np.vstack([xyz[:i], xyz[i+1:]])
                        d = calc_distance_pbc(xyz_i, xyz_j, box.length)
                        temp_g_r, _ = np.histogram(d, bins=n_bins, range=r_range)
                        g_r += temp_g_r

                # type_i-type_i
                elif pairs[0] == pairs[1]:
                    xyz_0 = xyz[types == pairs[0]]
                    for i, xyz_i in enumerate(xyz_0):
                        xyz_j = np.vstack([xyz_0[:i], xyz_0[i+1:]])
                        d = calc_distance_pbc(xyz_i, xyz_j, box.length)
                        temp_g_r, _ = np.histogram(d, bins=n_bins, range=r_range)
                        g_r += temp_g_r

                # type_i-type_j
                else:
                    for i, xyz_i in enumerate(xyz[types == pairs[0]]):
                        xyz_j = xyz[types == pairs[1]]
                        d = calc_distance_pbc(xyz_i, xyz_j, box.length)
                        temp_g_r, _ = np.histogram(d, bins=n_bins, range=r_range)
                        g_r += temp_g_r
            rho += (i + 1) / box.volume

    r = 0.5 * (edges[1:] + edges[:-1])
    V = 4./3. * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
    norm = rho * i
    g_r /= norm * V
    return r, g_r

