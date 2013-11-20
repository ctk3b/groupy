from groupy.gbb import *
from groupy.mdio import *
import pdb

sys = Gbb()
box = sys.load_lammps_data('example_inputs/data.peg6_0.2')
sys.wrap(box)
#sys.mirror(box)
write_xyz('example_outputs/peg6_0.2.xyz', sys.xyz, sys.types)
write_gro(sys, box,
    'example_outputs/peg6_0.2.gro', 'peg6_0.2_mirrored')
write_lammps_data(sys, box,
    'example_outputs/data.peg6_0.2_mirrored', 'peg6_0.2_mirrored')
