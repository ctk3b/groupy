#!/Users/CTK/canopy/bin/python
import sys

from groupy.gbb import *
from groupy.mdio import *

for name in sys.argv[1:]:
    sys = Gbb()
    sys.load_lammps_data('{0}'.format(name))

    write_xyz(sys.xyz, sys.types, '{0}.xyz'.format(name))

