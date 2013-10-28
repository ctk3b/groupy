from groupy.mdio import *
from groupy.order import *
from groupy.visualization import *

# xyz and mass data for a single MPC molecule
xyz, types = read_xyz('data_files/mpc.xyz')
mass = np.loadtxt('data_files/mpc_mass.txt')

I = calc_inertia_tensor(xyz, mass)
director = calc_director(I)

splat(xyz, director, types)

