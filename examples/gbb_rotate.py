from groupy.gbb import *
from groupy.mdio import *
from groupy.visualization import *


sys = Gbb()
box = sys.load_xyz('example_inputs/mpc.xyz')
splat(sys.xyz, sys.types)

sys.rotate([90.0, 0.0, 0.0])
splat(sys.xyz, sys.types)

sys.rotate([0.0, 90.0, 0.0])
splat(sys.xyz, sys.types)

sys.rotate([0.0, 0.0, 90.0])
splat(sys.xyz, sys.types)
