from groupy.gbb import *
from groupy.mdio import *
from groupy.visualization import *


sys = Gbb()
box = sys.load_xyz('example_inputs/mpc.xyz')
splat(sys.xyz, sys.types)

sys.translate([10.0, 10.0, 10.0])
splat(sys.xyz, sys.types)
