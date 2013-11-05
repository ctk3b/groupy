from groupy.gbb import *
from groupy.mdio import *
from groupy.visualization import *
from groupy.system import *


pegs = list()
for i in range(2):
    peg = Gbb()
    box = peg.load_lammps_data('example_inputs/data.peg6_0.2')
    pegs.append(peg)

pegs[0].translate(0, 0, 50)

sys = System()
sys.append_gbbs(pegs)

write_xyz('example_outputs/blah.xyz', sys.xyz, sys.types)
write_lammps_data(sys, box, 'example_outputs/data.blah')
