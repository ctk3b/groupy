from groupy.gbb import *
from groupy.mdio import *

sys = Gbb()
sys.load_data('data_files/data.peg6_0.2')

write_xyz('peg6_0.2.xyz', sys.xyz, sys.types)
