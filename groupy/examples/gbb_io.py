from groupy.gbb import *
from groupy.mdio import *

sys = Gbb()
box = sys.load_data('data_files/data.peg6_0.2')
sys.wrap(box)
sys.mirror(box)
write_xyz('peg6_0.2.xyz', sys.xyz, sys.types)
