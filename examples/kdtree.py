from groupy.gbb import *
from groupy.mdio import *
from groupy.visualization import *
from groupy.system import *
from groupy.general import *


'''centers = list()
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
'''
lt = 'unwrapped.lammptrj'

with open(lt, 'r') as trj:
	while True:
		try:
			xyz, types, step, box = read_frame_lammpstrj(trj)
		except:
			print 'End of file'
			break

		CER_coords = xyz[0,72*49]
		CER_molecules = np.split(CER_coords, range(0, 72*49, 49))
		CHOL_coords = 



neighbors = get_points_in_range(centers, centers[25], 200, 7)
splat(centers[neighbors], highlight=0, dims=dims)


