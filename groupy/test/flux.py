from groupy.general import *
from mpl_toolkits.mplot3d import Axes3D

# --- user input ---
coverage = 0.4
n_planes = 100  # number of flux planes
do_calc = True  # redo calculation, otherwise read from file
over_write_data = True

# generate intermediate plots
plot_vs_time = False
plot_vs_dist = False

# --- system info ---
if coverage == 0.2:
    # atom indices of regions
    botsub_ids = range(3200)
    botfilm_ids = range(3200, 4180)
    water_ids = range(4180, 6880) + range(11060, 13760)
    topsub_ids = range(6880, 10080)
    topfilm_ids = range(10080, 11060)
    # normal force vs distance data
    fn = np.array([7222.51282396, 5893.61521886, 4870.86290104, 3731.98751177, 2850.51780246,
        2173.26935712, 1595.29090095, 1181.89105584, 949.23856132, 586.38516427,
        399.97429248,  238.17032324, 134.06351658, 70.31856361])
    distances = range(25, 39)
elif coverage == 0.4:
    # atom indices of regions
    botsub_ids = range(3200)
    botfilm_ids = range(3200, 5060)
    water_ids = range(5060, 7160) + range(12220, 14320)
    topsub_ids = range(7160, 10360)
    topfilm_ids = range(10360, 12220)
    # normal force vs distance data
    fn = np.array([7631.08622339, 5984.56342842, 5007.35503764, 4378.8169841, 3469.1465996,
         2652.41144195, 1968.97817919, 1425.19978545, 877.62964466, 584.9258169,
         357.79081326,  232.29089716, 59.67694195, -44.29309985])
    #distances = range(28, 42)
    distances = range(28, 30)

system_info = {'botsub': botsub_ids,
               'botfilm': botfilm_ids,
               'water': water_ids,
               'topsub': topsub_ids,
               'topfilm': topfilm_ids}

fnvd = dict(zip(distances, fn))

# --- main ---
if do_calc:
    fn_z_flux = np.empty(shape=(len(distances), n_planes))
    all_planes = np.empty(shape=(len(distances), n_planes))
    for d, dist in enumerate(distances):
        file_name = ('/Users/CTK/Science/iMoDELS/data/mpc/shearing/'
            'mpc_' + str(coverage) + '_h2o_shear_' + str(dist) + '.lammpstrj')

        # read first frame to get system dimensions
        with open(file_name, 'r') as trj:
            xyz, _, _, _ = read_frame_lammpstrj(trj)

        # calculate surface area
        x_max = xyz[system_info['botsub']][:, 0].max()
        x_min = xyz[system_info['botsub']][:, 0].min()
        y_max = xyz[system_info['botsub']][:, 1].max()
        y_min = xyz[system_info['botsub']][:, 1].min()
        area = (x_max - x_min) * (y_max - y_min)

        # determine bounds on planes based on substrate surface coords
        z_min = xyz[system_info['botsub']][:, 2].max()
        z_max = xyz[system_info['topsub']][:, 2].min()
        planes = np.linspace(z_min, z_max, n_planes)

        fluxes, steps = calc_flux(file_name,
                system_info,
                planes,
                area,
                max_time=5e6)
        # convert atoms / (angstrom^2 * fs) to h2o_molecules / (nm^2 * ps)
        fluxes = fluxes * (100 * 1000 / 3)

        # show fluxes and averages over time
        if plot_vs_time:
            fig = plt.figure()
        # calculate time averaged flux across each plane
        avg_fluxes = np.empty(shape=(n_planes))
        for i in range(fluxes.shape[1]):
            p = np.polyfit(steps, fluxes[:, i], 1)
            avg_fluxes[i] = p[1]
            if plot_vs_time:
                plt.plot(steps, fluxes[:, i], '--', alpha=0.4)
                plt.plot([0, max(steps)], [p[1], p[1]])
        if plot_vs_time:
            plt.show()

        # store all the final data we want to keep
        all_planes[d, :] = planes
        fn_z_flux[d, :] = avg_fluxes

        # show fluxes at respective z-coords
        if plot_vs_dist:
            fig = plt.figure()
            plt.plot(avg_fluxes, planes)
            plt.xlabel(ur'Flux (molecules/nm$\mathregular{^2\cdot}$ ps)')
            plt.ylabel(ur'z (\u00c5)')
            fig.savefig('slice.pdf', bbox_inches='tight')

        print 'Done reading coverage %.1f at D=%2d ' % (coverage, dist)
    if over_write_data:
        np.savetxt('data_flux' + str(coverage) + '.txt', fn_z_flux)
        np.savetxt('data_planes' + str(coverage) + '.txt', all_planes)
else:
    fn_z_flux = np.loadtxt('data_flux' + str(coverage) + '.txt')
    all_planes = np.loadtxt('data_planes' + str(coverage) + '.txt')

# --- plotting ---
#plot_type = 'boring'
plot_type = 'imshow'
#plot_type = 'contourf'
#plot_type = '3D!!!'

# using imshow
if plot_type == 'imshow':
    z = np.arange(n_planes)

    fig, ax = plt.subplots()
    extent = [fn.min(), fn.max(), z.min(), z.max()]
    cax = ax.imshow(np.rot90(fn_z_flux, 3),
            extent=extent,
            origin='lower',
            aspect='auto')
    #for force in fn:
    #    ax.axvline(force, z.min(), z.max(), color='k', linestyle='--')
    #ax.axvline(1500.0, z.min(), z.max(),
    #        color='k',
    #        linestyle='--',
    #        linewidth=2)

    ax.set_xticks(range(0, 7001, 1000))
    ax.set_xlim([0, 7400])
    ax.set_xlabel(r'Normal force per area (MPa)')

    ax.set_yticklabels(np.arange(0.0, 1.1, 0.2))
    ax.set_ylim([z.min(), z.max() + 1])
    #ax.set_ylabel(r'Plane #')
    ax.set_ylabel(r'$\mathregular{\frac{z\ -\ z_{min}}{z_{max}\ -\ z_{min}}}$',
        fontsize=20)

    cb = fig.colorbar(cax, ticks=np.arange(0, 0.17, 0.02))
    cb.set_label(r'Flux (molecules/[nm$\mathregular{^2\cdot}$ ps])')

    fig.savefig('test_imshow.pdf', bbox_inches='tight')

# using contourf
if plot_type == 'contourf':
    z = np.arange(n_planes)

    fig, ax = plt.subplots()
    X, Y = np.meshgrid(fn, z)
    Z = fn_z_flux.T

    cax = ax.contourf(X, Y, Z)

    ax.set_xticks(range(0, 7001, 1000))
    ax.set_xlim([0, 7400])
    ax.set_xlabel(r'Normal force per area (MPa)')

    ax.set_yticklabels(np.arange(0.0, 1.1, 0.2))
    ax.set_ylim([z.min(), z.max() + 1])
    ax.set_ylabel(r'$\mathregular{\frac{z\ -\ z_{min}}{z_{max}\ -\ z_{min}}}$',
        fontsize=20)

    cb = fig.colorbar(cax, ticks=np.arange(0, 0.17, 0.02))
    cb.set_label(r'Flux (molecules/[nm$\mathregular{^2\cdot}$ ps])')

    fig.savefig('test_contour.pdf', bbox_inches='tight')

# using plane ole lines
if plot_type == 'boring':
    # re-formulate z-values as distance to wall
    d_from_walls = np.copy(all_planes[::2])
    for i, planes in enumerate(all_planes[::2]):
        mid = (max(planes) + min(planes)) / 2
        d_from_walls[i] = np.array(planes) - mid

    fig, ax = plt.subplots()
    for i, force in enumerate(fn[::2]):
        ax.plot(d_from_walls[i], fn_z_flux[::2][i],
                label='$\mathregular{F_n}$=%.0f' % (force))
    ax.legend(bbox_to_anchor=(1.3, 1.0))

    fig.savefig('test_boring.pdf', bbox_inches='tight')

# using 3D projections
if plot_type == '3D!!!':
    # re-formulate z-values as distance to wall
    d_from_walls = np.copy(all_planes)
    for i, planes in enumerate(all_planes):
        mid = (max(planes) + min(planes)) / 2
        d_from_walls[i] = np.array(planes) - mid

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for i, force in enumerate(fn):
        v_force = np.empty(d_from_walls[i].shape[0])
        v_force.fill(force)
        ax.scatter(v_force, d_from_walls[i], zs=fn_z_flux[i],
                label='$\mathregular{F_n}$=%.0f' % (force))

    ax.set_xlabel(r'Normal force per area (MPa)')

    ax.set_ylabel(r'Distance from center ($\mathregular{\AA}$)')

    ax.set_zlabel(r'Flux (molecules/[nm$\mathregular{^2\cdot}$ ps])')


    fig.savefig('test_3D.pdf', bbox_inches='tight')

