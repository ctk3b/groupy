import pdb

from groupy.gbb import *
from groupy.mdio import *
from groupy.visualization import *


sys = Gbb()
box = sys.load_xyz('example_inputs/mpc.xyz')
splat(sys.xyz, sys.types)
sys.rotate([True, False, False], [90.0, 0.0, 0.0])
splat(sys.xyz, sys.types)
sys.rotate([False, True, False], [0.0, 90.0, 0.0])
splat(sys.xyz, sys.types)
sys.rotate([False, False, True], [0.0, 0.0, 90.0])
splat(sys.xyz, sys.types)
write_xyz('example_outputs/mpc_rotated.xyz', sys.xyz, sys.types)
