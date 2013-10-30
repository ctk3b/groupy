import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pdb

# type: (color, vdw radius)
a_info = {'C': ('teal', 1.7),
          'H': ('white', 1.2),
          'O': ('red', 1.52),
          'N': ('blue', 1.55),
          'P': ('orange', 1.8),
          'Si': ('yellow', 2.1)}
a_info_maya = {'C': ((0, 0.8, 0.8), 1.7),
               'H': ((1, 1, 1), 1.2),
               'O': ((1, 0, 0), 1.52),
               'N': ((0, 0, 1), 1.55),
               'P': ((1, 0.5, 0), 1.8),
               'Si': ((1, 1, 0), 2.1)}
a_info_maya = {'C': (2, 1.7),
               'H': (1, 1.2),
               'O': (3, 1.52),
               'N': (4, 1.55),
               'P': (5, 1.8),
               'Si': (6, 2.1)}




def splat(xyz, types=None, direction=None):
    """Dump coordinates into 3D plot and show it using matplotlib.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i, atom in enumerate(xyz):
        if types is None:
            ax.scatter(atom[0], atom[1], atom[2],
            marker='o',
            s=100)
        else:
            ax.scatter(atom[0], atom[1], atom[2],
            c=a_info[types[i]][0],
            marker='o',
            s=100 * a_info[types[i]][1] ** 3)


    bounds = [xyz.min(), xyz.max()]
    ax.set_xlim(bounds)
    ax.set_ylim(bounds)
    ax.set_zlim(bounds)

    if direction is None:
        pass
    else:
        # TODO: robust way to choose origin
        ax.plot([0, 20*direction[0]], [0, 20*direction[1]], [0, 20*direction[2]], 
            c='k', linewidth=10)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

def splat_maya(xyz, types=None, direction=None):
    """Dump coordinates into 3D plot and show it using mayavi.

    TODO:
        -resolve crashing window issue
        -choosing specific colors instead of by scalar
        -show director
    """
    import mayavi.mlab as mlab
    if types is None:
        mlab.points3d(xyz[:, 0], xyz[:, 1], xyz[:, 2])
    else:
        colors = [a_info_maya[x][0] for x in types]
        sizes = [a_info_maya[x][1] for x in types]
        pts = mlab.quiver3d(xyz[:, 0], xyz[:, 1], xyz[:, 2],
                sizes, sizes, sizes, scalars=colors, mode='sphere')
        pts.glyph.color_mode = 'color_by_scalar'
        pts.glyph.glyph_source.glyph_source.center = [0, 0, 0]
    
    mlab.show()
