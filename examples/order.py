from groupy.mdio import *
from groupy.box import Box
from groupy.system import System
from groupy.order import *
import pdb

info = [(36, 74, 'chol')]
counter = 0
atoms = [i for i in range(24)]
pdb.set_trace()

s2 = []
with open('example_inputs/monolayer600K.lammpstrj', 'r') as trj:
    while True and counter < 10:
        try:
            xyz, types, step, box = read_frame_lammpstrj(trj)
            system = System(system_info=info, box=box)
            system.convert_from_traj(xyz, types)
        except:
            print "Reached end of file"
            break

        directors = []
        for lipid in system.gbbs:
            lipid.masses = np.ones((lipid.xyz.shape[0]))
            #print lipid.calc_inertia_tensor()
            directors.append(calc_director(lipid.calc_inertia_tensor(atoms)))
        s2.append(calc_S2(calc_Q_tensor(np.asarray(directors))))
        counter += 1

print np.mean(s2)
