from groupy.gbb import *
from groupy.mdio import *
import pdb

sys = Gbb()
dims = sys.load_lammps_data('example_inputs/data.peg6_0.2')

box = Box(dims[:, 1], dims[:, 0])
sys.wrap(box)
#sys.mirror(box)
write_xyz(sys.xyz, sys.types, 'example_outputs/peg6_0.2.xyz')
write_lammps_data(sys, box,
    'example_outputs/data.peg6_0.2_mirrored', 'peg6_0.2_mirrored')
