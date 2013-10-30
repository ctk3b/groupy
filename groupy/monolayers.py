import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.stats
import scipy.integrate

from groupy.mdio import *


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
        heights = np.empty(shape=(2, 0))
        while True:
            try:
                xyz, _, step, _ = read_frame_lammpstrj(trj)
            except:
                print "Reached end of '" + file_name + "'"
                break
            # temp container for z-coords of top, [0], and bottom, [1], films
            temp_heights = np.empty(shape=(2, len(system_info['botfilm'])))
            temp_heights[0] = xyz[system_info['topfilm']][:, 2]
            temp_heights[1] = xyz[system_info['botfilm']][:, 2]
            heights = np.hstack((heights, temp_heights))

        top_bounds = find_cutoff('top', heights, plot=True)
        bot_bounds = find_cutoff('bot', heights, plot=True)
    return top_bounds, bot_bounds


def calc_flux(file_name,
        system_info,
        planes,
        area,
        max_time=np.inf):
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


def calc_res_time(file_name,
        system_info,
        top_bounds,
        bot_bounds,
        slab,
        max_time=np.inf,
        return_data=False,
        plot=False):
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
        ax1.hist(heights[i], color='green')
        ax1.vlines(z[idx], 0, 350, linewidth=3)

        ax2 = ax1.twinx()
        ax2.plot(z, cdf, color='blue', linewidth=3)
        pdf = sp.stats.rayleigh.pdf(z, loc=param[0], scale=param[1])
        ax2.plot(z, pdf, color='red', linewidth=3)

        fig.savefig(film + '_film_thickness.pdf', bbox_inches='tight')
    return film_bounds
