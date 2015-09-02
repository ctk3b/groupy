from collections import defaultdict
import pdb

import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace
import scipy as sp
import scipy.stats
import scipy.integrate
from scitools.numpyutils import meshgrid

from groupy.mdio import read_frame_lammpstrj
from groupy.general import find_nearest


def calc_vel_profile(file_name, system_info):
    """
    """
    with open(file_name, 'r') as trj:
        z_vx = dict.fromkeys(system_info.keys(), np.empty(shape=(0,2)))
        step = -np.Inf
        while True:
            try:
                xyz, _, step, _ = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '" + file_name + "'"
                break
            if step % 10000 == 0:
                print "Read step " + str(step)
            if step > 0:
                for region, indices in system_info.iteritems():
                    coords = xyz[indices]
                    prev_coords = prev_xyz[indices]
                    temp = np.zeros(shape=(coords.shape[0], 2))
                    # average z
                    temp[:, 0] = 0.5 * (coords[:, 2] + prev_coords[:, 2])
                    # x-velocity
                    temp[:, 1] = ((coords[:, 0] - prev_coords[:, 0])
                            / (step - prev_step))
                    z_vx[region] = np.vstack((z_vx[region], temp))

            prev_xyz = xyz
            prev_step = step
    return z_vx


def calc_film_heights(file_name, system_info):
    """Calculate z-coordinate bounds of top and bottom monolayers.

    Bounds are the atom closest to the substrate and the z-coordinate
    that on average includes 90% of the atoms throughout the simulation.

    Args:
        file_name (str): name of trajectory to read
        system_info (dict): dictionary containing indices of atoms in top and
            bottom monolayer. Corresponding keys must be:
            'topfilm' and 'botfilm'
    Returns:
        top_bounds (tuple): z-bounds of top monolayer (min, max)
        bot_bounds (tuple): z-bounds of bot monolayer (min, max)
    """
    with open(file_name, 'r') as trj:
        top_bounds = list()
        bot_bounds = list()
        while True:
            try:
                xyz, _, step, _ = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '" + file_name + "'"
                break

            # temp container for z-coords of top, [0], and bottom, [1], films
            heights = np.empty(shape=(2, len(system_info['botfilm'])))
            heights[0] = xyz[system_info['topfilm']][:, 2]
            heights[1] = xyz[system_info['botfilm']][:, 2]

            top = find_cutoff('top', heights)
            bot = find_cutoff('bot', heights, plot=True)
            top_bounds.append(top)
            bot_bounds.append(bot)
    return np.asarray(top_bounds), np.asarray(bot_bounds)


def calc_flux(file_name, system_info, planes, area, max_time=np.inf):
    """Calculate fluxes of water across multiple x-y planes.

    Args:
        file_name (str): name of trajectory to read
        system_info (dict): dictionary containing indices of water atoms
            Corresponding key must be: 'water'
        planes (np.ndarray): z-coords of planes to calculate flux through
        area (float): surface area of planes
    Returns:
        fluxes_over_time (np.ndarray): fluxes through 'planes'
    """
    steps = list()
    fluxes_over_time = np.empty(shape=(0, len(planes)))
    with open(file_name, 'r') as trj:
        step = -1
        while step < max_time:
            try:
                xyz, _, step, _ = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '" + file_name + "'"
                break

            # select z-coords of water atoms
            water = xyz[system_info['water']][:, 2]
            if step > 0:
                steps.append(step)
                fluxes = np.empty(shape=(1, len(planes)))
                for i, plane in enumerate(planes): # TODO: vectorize
                    # select water atoms that were and are above the flux plane
                    were_above = np.where(prev_water > plane)[0]
                    are_above = np.where(water > plane)[0]
                    # count how many left the level
                    n_fluxed = (len(were_above) -
                               len(np.intersect1d(are_above, were_above)))
                    # calc dat flux
                    fluxes[0, i] = n_fluxed / (area * (step - prev_step))
                fluxes_over_time = np.vstack((fluxes_over_time, fluxes))
            # store current frame
            prev_water = water
            prev_step = step
    return fluxes_over_time, steps


def calc_density(file_name, system_info, group, masses, axis, planes,
        area, max_frames=np.inf):
    """Calculate density along a specified axis.

    Args:
        file_name (str): name of trajectory to read
        system_info (dict): dictionary containing indices of groups in system
        group (list): list of keys in system_info
        masses (dict):  {atom_type: mass}
        axis (int): axis along which to calculate density
        planes (np.ndarray): coords of planes to calculate flux through
        area (float): area of plane normal to specified axis
        max_frames (int): maximum number of frames to read
    Returns:
        densities_over_time (np.ndarray): densities along axis over time
    """
    if max_frames == np.inf:
        densities_over_time = np.empty(len(planes) - 1)
    else:
        densities_over_time = np.empty(shape=(max_frames, len(planes) - 1))

    with open(file_name, 'r') as trj:
        step = -1
        while step < max_frames - 1:
            try:
                xyz, types, time_step, box = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '{0}'".format(file_name)
                break
            step += 1
            print "Read frame #{0}".format(step)
            planes[0] = box.mins[0]
            planes[-1] = box.maxs[0]

            # select coords of relevant atoms along specified axis
            # TODO: multiple groups
            selected = xyz[system_info[group]][:, axis]
            selected_types = types[system_info[group]]

            volumes = area * (np.diff(planes))

            if max_frames == np.inf:
                temp_masses = np.empty(len(planes) - 1)
            for i, plane in enumerate(planes[:-1]):
                # select atoms in layer
                atoms_in_layer = np.where((selected > plane)
                                     & (selected < planes[i+1]))[0]
                types_in_layer = selected_types[atoms_in_layer]
                # TODO: REMOVE HARDCODED TYPE OFFSET
                mass_in_layer = np.sum(masses[types_in_layer - 20])
                if max_frames == np.inf:
                    temp_masses[i] = mass_in_layer / volumes[i]
                else:
                    densities_over_time[step, i] = mass_in_layer / volumes[i]
            if max_frames == np.inf:
                densities_over_time = np.hstack((densities_over_time, temp_masses))
    return densities_over_time

def slab_density(file_name, system_info, group, masses, type_offset, axis,
        planes, area, max_frames=np.inf):
    """Calculate density in a slab.
    """
    with open(file_name, 'r') as trj:
        step = -1
        while step < max_frames - 1:
            try:
                xyz, types, time_step, box = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '{0}'".format(file_name)
                break
            step += 1
            print "Read frame #{0}".format(step)
            planes[0] = box.mins[0]
            planes[-1] = box.maxs[0]

            # select coords of relevant atoms along specified axis
            selected = xyz[system_info[group]][:, axis]
            selected_types = types[system_info[group]]

            volumes = area * (np.diff(planes))

            if max_frames == np.inf:
                temp_masses = np.empty(len(planes) - 1)
            for i, plane in enumerate(planes[:-1]):
                # select atoms in layer
                atoms_in_layer = np.where((selected > plane)
                                     & (selected < planes[i+1]))[0]
                types_in_layer = selected_types[atoms_in_layer]
                mass_in_layer = np.sum(masses[types_in_layer - type_offset])
                if max_frames == np.inf:
                    temp_masses[i] = mass_in_layer / volumes[i]
                else:
                    densities_over_time[step, i] = mass_in_layer / volumes[i]
            if max_frames == np.inf:
                densities_over_time = np.hstack((densities_over_time, temp_masses))
    return densities_over_time


def calc_res_time(file_name, system_info, top_bounds, bot_bounds, slab,
        max_time=np.inf, return_data=False, plot=False):
    """Calculate residence time of water molecules in monolayers.

    Args:
        file_name (str): name of trajectory to read
        system_info (dict): dictionary containing indices of water atoms
            Corresponding key must be: 'water'
        top_bounds (tuple): z-bounds of top monolayer (min, max)
        bot_bounds (tuple): z-bounds of bot monolayer (min, max)
        slab (tuple): bounds of slab containing initial water molecules
        return_data (bool): return time and R(t) values
        plot (bool): optional flag to plot fitted exponential decay
    Returns:
        res_time (float): residence time
        time (list): simulation times of frames read
        frac_remaining (list): fraction of water remaining in monolayer

    TODO:
        -multiple slabs simultaneously
    """
    with open(file_name, 'r') as trj:
        data = list()
        steps = list()

        # monolayer cutoff
        top_bound = top_bounds[0]
        bot_bound = bot_bounds[1]

        # define bounds of slab containing starting positions
        top_top_slab = top_bounds[1] - slab[0]
        top_bot_slab = top_bounds[1] - slab[1]
        bot_top_slab = bot_bounds[0] + slab[1]
        bot_bot_slab = bot_bounds[0] + slab[0]

        step = -1
        while step < max_time:
            try:
                xyz, _, step, _ = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '" + file_name + "'"
                break
            steps.append(step)
            if step == 0:
                # select z-coords of water atoms
                first = xyz[system_info['water']][:, 2]
                # select water atoms initially in top and bottom slabs
                start_in_topslab = np.where(((first > top_bot_slab)
                        & (first < top_top_slab)))[0]
                start_in_botslab = np.where(((first < bot_top_slab)
                        & (first > bot_bot_slab)))[0]
                # count em
                n_init = float(len(start_in_topslab)
                        + len(start_in_botslab))
                data.append(n_init)
            else:
                current = xyz[system_info['water']][:, 2]
                # select water atoms initially in slabs and still in respective monolayers
                still_in_toplayer = start_in_topslab[current[start_in_topslab] > top_bound]
                still_in_botlayer = start_in_botslab[current[start_in_botslab] < bot_bound]
                # count em
                diff = len(np.intersect1d(start_in_topslab, still_in_toplayer))
                diff += len(np.intersect1d(start_in_botslab, still_in_botlayer))
                data.append(diff)
    # convert to ps
    time = np.array([x / 1000. for x in steps])
    # normalize by number of atoms
    frac_remaining = np.array([x / n_init for x in data])

    # non-linear fit
    A0 = frac_remaining.max()
    K0 = -y.max() / t.max()
    C0 = np.mean(y[int(0.75 * len(y)):])

    guesses = [A0, K0, C0]
    A, K, C = fit_exp_nonlinear(time, frac_remaining, guesses)
    res_time = K

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        fit_y = model_exp(time, A, K, C)
        plot_fit(ax, time, frac_remaining, fit_y)
        fig.savefig('res_time_exp_fit.pdf', bbox_inches='tight')

    #res_time = sp.integrate.simps(frac_remaining, time)  # integrate RTCF

    if return_data:
        return res_time, time, frac_remaining
    else:
        return res_time


def find_cutoff(film, heights, plot=False):
    """Find z-coordinates that includes 90% of the points in 'heights'.

    Fits a Rayleigh distribution function a row 'heights' and calculates the
    z-coordinates that encompass 90% of the cumulative distribution function.

    NOTE: 'top' z_max is the atom closest to the top substrate.
          'bot' z_max is the atom farthest from the bottom substrate.
    This simply means we're interested in the top 90% of the distribution
    distribtution when considering the top film and bottom 90% when considering
    the bottom film.

    Args:
        film (str): film identifier. Must be 'top' or 'bot'
        heights (numpy.ndarray): array containing z-coords of the top, [0],
                                 and bottom, [1], film atoms
        plot (bool): optional flag to print histogram and fitted distributions

    Returns:
        film_bounds (tuple): bounds of monolayer film
    """
    if film == 'top':
        i = 0  # indexing of 'heights'
        cut = 0.1
    elif film == 'bot':
        i = 1
        cut = 0.9
    else:
        raise Exception("Invalid film identifier! Must be 'top' or 'bot'.")
    z_min = heights[i].min()
    z_max = heights[i].max()

    # fit Rayleigh distribution function
    z = np.linspace(z_min, z_max, 200)
    param = sp.stats.rayleigh.fit(heights[i])
    cdf = sp.stats.rayleigh.cdf(z, loc=param[0], scale=param[1])

    # find 90% cutoff value
    idx, val = find_nearest(cdf, cut)
    if film == 'top':
        film_bounds = (z[idx], z_max)
    elif film == 'bot':
        film_bounds = (z_min, z[idx])

    if plot:
        fig, ax1 = plt.subplots()
        ax1.set_xlabel(u'z (\u00c5)')
        ax1.set_ylabel('Count')
        ax1.hist(heights[i], color='green')
        ax1.vlines(z[idx], 0, 350, linewidth=3)

        ax2 = ax1.twinx()
        ax2.set_ylabel('Cumulative distribution function')
        ax2.plot(z, cdf, color='blue', linewidth=3)
        pdf = sp.stats.rayleigh.pdf(z, loc=param[0], scale=param[1])
        ax2.plot(z, pdf, color='red', linewidth=3)

        fig.savefig(film + '_film_thickness.pdf', bbox_inches='tight')
    return film_bounds


def voxel_density(file_name, system_info, box, n_grid=[50, 50, 10], z_bounds=[], max_time=np.Inf):
    """
    """

    count = defaultdict(int)
    n_frames = 0
    x_min = box.mins[0]
    x_max = box.maxs[0]
    y_min = box.mins[1]
    y_max = box.maxs[1]
    z_min = z_bounds[0]
    z_max = z_bounds[1]

    xs = linspace(x_min, x_max, n_grid[0])
    ys = linspace(y_min, y_max, n_grid[1])
    zs = linspace(z_min, z_max, n_grid[2])

    vol = (x_max-x_min) * (y_max-y_min) * (z_max-z_min)
    vol_per_voxel = vol / np.prod(n_grid)
    print 'Volume of voxel: {0}'.format(vol_per_voxel)

    units = 1.660538  # au/ang^3 to g/cm^3
    mass_per_volume = {1: 1.008 / vol_per_voxel * units,
            2: 14.007 / vol_per_voxel * units,
            3: 12.011 / vol_per_voxel * units,
            4: 12.011 / vol_per_voxel * units,
            5: 1.008 / vol_per_voxel * units,
            6: 12.011 / vol_per_voxel * units,
            7: 1.008 / vol_per_voxel * units,
            8: 15.999 / vol_per_voxel * units,
            9: 30.974 / vol_per_voxel * units,
            10: 15.990 / vol_per_voxel * units,
            11: 12.011 / vol_per_voxel * units,
            12: 15.999 / vol_per_voxel * units,
            13: 12.011 / vol_per_voxel * units,
            14: 15.999 / vol_per_voxel * units,
            15: 12.011 / vol_per_voxel * units,
            16: 12.011 / vol_per_voxel * units,
            17: 1.008 / vol_per_voxel * units,
            18: 12.011 / vol_per_voxel * units,
            19: 15.999 / vol_per_voxel * units,
            20: 1.008 / vol_per_voxel * units,
            21: 28.085 / vol_per_voxel * units,
            22: 28.085 / vol_per_voxel * units,
            23: 15.999 / vol_per_voxel * units,
            24: 1.008 / vol_per_voxel * units,
            25: 28.085 / vol_per_voxel * units,
            26: 15.999 / vol_per_voxel * units,
            27: 15.999 / vol_per_voxel * units,
            28: 1.008 / vol_per_voxel * units}

    with open(file_name, 'r') as trj:
        step = -np.Inf
        while step < max_time:
            try:
                xyz, types, _, box = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '" + file_name + "'"
                break
            n_frames += 1

            #all_water = xyz[system_info['water']]
            #water_ids = np.where((all_water[:, 2] > z_min) & (all_water[:, 2] < z_max))
            #water = all_water[water_ids]
            bounded_ids = np.where((xyz[:, 2] > z_min) & (xyz[:, 2] < z_max))
            bounded_atoms = xyz[bounded_ids]

            #all_water_types = types[system_info['water']]
            #water_types = all_water_types[water_ids]
            bounded_types = types[bounded_ids]

            for atom, a_type in zip(bounded_atoms, bounded_types):
                # wrap coord if necessary
                for k, c in enumerate(atom):
                    if c < box.dims[k, 0]:
                        atom[k] = box.dims[k, 1] - abs(box.dims[k, 0] - c)
                    elif c > box.dims[k, 1]:
                        atom[k] = box.dims[k, 0] + abs(c - box.dims[k, 1])
                x = np.where(np.histogram([atom[0]], bins=xs)[0] ==  1)[0][0]
                y = np.where(np.histogram([atom[1]], bins=ys)[0] ==  1)[0][0]
                z = np.where(np.histogram([atom[2]], bins=zs)[0] ==  1)[0][0]

                count[(x, y, z)] += mass_per_volume[a_type]

    for voxel, density in count.items():
        count[voxel] = density / n_frames
    return count, vol_per_voxel


def pore_distribution(file_name,
        system_info,
        groups,
        bounds,
        n_bins=100,
        max_time=np.inf):
    """Calculate distribution of species across 2D pore

    Args:
        file_name (str): name of trajectory to read
        system_info (dict): dictionary containing indices of water atoms
            Corresponding key must be: 'water'
        groups (dict): dictionary of types that should be treated as a group
    Returns:
        densities_over_time (np.ndarray): densities along z-axis over time
    """

    counts = dict()
    for group in groups:
        counts[group], edges = np.histogram([0], bins=n_bins, range=(bounds[0], bounds[1]))
        counts[group][0] = 0

    print "Reading '" + file_name + "'"
    with open(file_name, 'r') as trj:
        step = -1
        while step < max_time:
            try:
                xyz, _, step, _ = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '" + file_name + "'"
                break

            for group in groups:
                temp_xyz = xyz[system_info[group]][:, 2]
                temp, _ = np.histogram(temp_xyz, bins=n_bins, range=(bounds[0], bounds[1]))
                counts[group] += temp

    return counts, edges

