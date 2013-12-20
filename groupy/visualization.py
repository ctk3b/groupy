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
"""
# corresponding RGB values for colors
a_info_maya = {'C': ((0, 0.8, 0.8), 1.7),
               'H': ((1, 1, 1), 1.2),
               'O': ((1, 0, 0), 1.52),
               'N': ((0, 0, 1), 1.55),
               'P': ((1, 0.5, 0), 1.8),
               'Si': ((1, 1, 0), 2.1)}
"""
# simplified color values to work with color_by_scalar
a_info_maya = {'C': (2, 1.7),
               'H': (1, 1.2),
               'O': (3, 1.52),
               'N': (4, 1.55),
               'P': (5, 1.8),
               'Si': (6, 2.1)}


def splat(xyz, types=None, direction=None, highlight=None, dims=None):
    """Dump coordinates into 3D plot and show it using matplotlib.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if types is None:
        ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2],
                marker='o',
                s=100)
    else:
        colors = [a_info[x][0] for x in types]
        sizes = [a_info[x][1] for x in types]
        ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2],
                c=colors,
                marker='o',
                s=[10 * x ** 3 for x in sizes])

    if highlight is not None:
        ax.scatter(xyz[highlight][0], xyz[highlight][1], xyz[highlight][2],
                marker='D',
                c='red',
                s=200)

    if dims is None:
        bounds = [xyz.min(), xyz.max()]
        ax.set_xlim(bounds)
        ax.set_ylim(bounds)
        ax.set_zlim(bounds)
    else:
        assert dims.shape == (3, 2)
        ax.set_xlim(dims[0, 0], dims[0, 1])
        ax.set_ylim(dims[1, 0], dims[1, 1])
        ax.set_zlim(dims[2, 0], dims[2, 1])

    if direction is None:
        pass
    else:
        assert direction.shape[0] == 3
        # TODO: robust way to choose origin
        zl = xyz[:, 2].argmin()
        zh = xyz[:, 2].argmax()
        d = np.linalg.norm(zh - zl)
        '''
        ax.plot([xyz[zl, 0], d * direction[0]],
                [xyz[zl, 1], d * direction[1]],
                [xyz[zl, 2], d * direction[2]],
                c='k', linewidth=5)
        '''
        ax.plot([0, d * direction[0]],
                [0, d * direction[1]],
                [0, d * direction[2]],
                c='k', linewidth=5)



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
    try:
        import mayavi.mlab as mlab
    except:
        print ('Failed to import mayavi. '
            + 'Please ensure that mayavi is properly installed.')
        return

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
