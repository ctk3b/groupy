from groupy.mdio import *
from groupy.inertia import *
from groupy.visualization import *

xyz, types = read_xyz('mpc.xyz')
mass = np.loadtxt('mpc_mass.txt')

I = calc_inertia_tensor(xyz, mass)
director = calc_vector(I)

splat(xyz, director, types)

