from groupy.order import *
from groupy.visualization import *
from groupy.gbb import *

mpc = Gbb()
mpc.load_xyz('example_inputs/mpc.xyz')
mpc.load_mass('example_inputs/mpc_mass.txt')

I = mpc.calc_inertia_tensor()
director = calc_director(I)

# three different splat modes
splat(mpc.xyz, mpc.types, director)  # with the director
splat(mpc.xyz, mpc.types)  # without the director
splat(mpc.xyz)  # without pretty colors and sizing
