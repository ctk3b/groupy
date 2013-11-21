import pdb

from groupy.gbb import *
from groupy.mdio import *
from groupy.visualization import *
from groupy.system import *
from groupy.general import *


centers = list()
mpcs = list()
# create a bunch of chains in space
for i in range(4):
    for j in range(4):
        for k in range(4):
            mpc = Gbb()
            box = mpc.load_xyz('example_inputs/mpc.xyz')
            mpc.load_mass('example_inputs/mpc_mass.txt')
            mpc.translate(i*30, j*30, k*30)
            # calculate and store all of the centers of mass
            mpc.calc_com()
            centers.append(mpc.com)
            mpcs.append(mpc)

dims=np.array([[0, 100], [0, 100], [0, 100]])

# show the chains
sys = System()
sys.append_gbbs(mpcs)
splat(sys.xyz, sys.types, dims=dims)

# show the centers of mass and highlight the one we want to look at
centers = np.array(centers)
splat(centers, highlight=25, dims=dims)

# show 6 closest neighbors to com of chain # 25
neighbors = get_points_in_range(centers, centers[25], 200, 7)
splat(centers[neighbors], highlight=0, dims=dims)


