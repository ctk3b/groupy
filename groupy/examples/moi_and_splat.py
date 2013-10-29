import pdb

from groupy.mdio import *
from groupy.order import *
from groupy.visualization import *
from groupy.gbb import *

mpc = gbb()
mpc.load_xyz('data_files/mpc.xyz')
mpc.load_mass('data_files/mpc_mass.txt')

I = mpc.calc_inertia_tensor()
director = calc_director(I)

# three different splat modes
splat(mpc.xyz, mpc.types, director)
splat(mpc.xyz, mpc.types)
splat(mpc.xyz)

