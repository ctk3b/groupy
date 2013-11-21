from groupy.gbb import *
from groupy.system import *
from groupy.mdio import *
from groupy.visualization import *
from groupy.box import *

def can_add(gbb, gbb_list, box, r_cut):
    """
    """
    for test_atom in gbb.xyz:
        for item in gbb_list:
            r_ij = calc_distance_pbc(item.xyz, test_atom, box.maxs)
            if (r_ij < r_cut).any():
                return False
    return True

peg = Gbb()
peg.load_mass("example_inputs/peg6_mass.txt")
peg.load_xyz("example_inputs/peg6.xyz")
peg.shift_com_to_origin()

box = Box(100)
pegs = list()

for i in range(100):
    t_peg = copy.copy(peg)
    t_peg.rotate(180 * np.random.rand(3))
    t_peg.translate(box.length * np.random.rand(3))
    if can_add(t_peg, pegs, box, 3.0):
        pegs.append(t_peg)
        print "Added " + str(i)

sys = System(box)
sys.append_gbbs(pegs)
splat(sys.xyz, sys.types, dims=box.dims)
write_xyz('example_outputs/random_pegs.xyz', sys.xyz, sys.types)
