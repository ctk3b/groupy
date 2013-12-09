import matplotlib.pyplot as plt
from scipy.integrate import simps

from groupy.rdf import *

traj_file = 'example_inputs/water.lammpstrj'
r, g_r = calc_rdf(traj_file, n_bins=80, max_frames=1)
#r, g_r = calc_rdf(traj_file, pair_types=[1, 1], n_bins=80)
#r, g_r = calc_rdf(traj_file, pair_types=[1, 2], n_bins=80)

#r, g_r = calc_rdf(traj_file, pair_types=[(1, 1)], n_bins=80, opencl=True)

plt.plot(r, g_r, 'ro-')
plt.plot([0, 8], [1, 1], 'k-')
plt.xlabel(u'r (\u212B)')
plt.ylabel('g(r)')
plt.show()
