import numpy as np


def get_2d_block(dx=0.01, length=1.0, height=1.0, center=np.array([0., 0.])):
    n1 = int(length / dx) + 1
    n2 = int(height / dx) + 1
    x, y = np.mgrid[-length / 2.:length / 2.:n1 *
                    1j, -height / 2.:height / 2.:n2 * 1j]
    x, y = np.ravel(x), np.ravel(y)
    return x + center[0], y + center[1]


def hydrostatic_tank_2d(fluid_length=1., fluid_height=2.,
                        tank_height=2.3, tank_layers=2,
                        fluid_spacing=0.1, tank_spacing=0.1, close=False):
    xf, yf = get_2d_block(dx=fluid_spacing,
                          length=fluid_length,
                          height=fluid_height)

    xt_4 = np.array([])
    yt_4 = np.array([])
    if close == False:
        xt_1, yt_1 = get_2d_block(dx=fluid_spacing,
                                  length=tank_layers*fluid_spacing,
                                  height=tank_height+fluid_spacing/2.)
        xt_1 -= max(xt_1) - min(xf) + fluid_spacing
        yt_1 += min(yf) - min(yt_1)

        xt_2, yt_2 = get_2d_block(dx=fluid_spacing,
                                  length=tank_layers*fluid_spacing,
                                  height=tank_height+fluid_spacing/2.)
        xt_2 += max(xf) - min(xt_2) + fluid_spacing
        yt_2 += min(yf) - min(yt_2)

    else:
        xt_1, yt_1 = get_2d_block(dx=fluid_spacing,
                                  length=tank_layers*fluid_spacing,
                                  height=fluid_height + tank_layers * fluid_spacing)
        xt_1 -= max(xt_1) - min(xf) + fluid_spacing
        yt_1 += min(yf) - min(yt_1)

        xt_2, yt_2 = get_2d_block(dx=fluid_spacing,
                                  length=tank_layers*fluid_spacing,
                                  height=fluid_height + tank_layers * fluid_spacing)
        xt_2 += max(xf) - min(xt_2) + fluid_spacing
        yt_2 += min(yf) - min(yt_2)

        xt_3, yt_3 = get_2d_block(dx=fluid_spacing,
                                length=max(xt_2) - min(xt_1),
                                height=tank_layers*fluid_spacing)
        yt_3[:] = yt_3[:] - (max(yt_3) - min(yf)) - fluid_spacing

        xt_4, yt_4 = get_2d_block(dx=fluid_spacing,
                                length=max(xt_2) - min(xt_1),
                                height=tank_layers*fluid_spacing)
        yt_4[:] = yt_4[:] + max(yf) - min(yt_4) + fluid_spacing

    xt_3, yt_3 = get_2d_block(dx=fluid_spacing,
                              length=max(xt_2) - min(xt_1),
                              height=tank_layers*fluid_spacing)
    yt_3[:] = yt_3[:] - (max(yt_3) - min(yf)) - fluid_spacing

    xt = np.concatenate([xt_1, xt_2, xt_3, xt_4])
    yt = np.concatenate([yt_1, yt_2, yt_3, yt_4])

    # plt.scatter(xt_3, yt_3)
    # plt.show()

    return xf, yf, xt, yt


def translate_system_with_left_corner_as_origin(x, y, z):
    translation = [min(x), min(y), min(z)]
    x[:] = x[:] - min(x)
    y[:] = y[:] - min(y)
    z[:] = z[:] - min(z)
    return translation


def create_tank_2d_from_block_2d(xf, yf, tank_length, tank_height,
                                 tank_spacing, tank_layers, close=False):
    """
    This is mainly used by granular flows

    Tank particles radius is spacing / 2.
    """
    ####################################
    # create the left wall of the tank #
    ####################################
    xleft, yleft = get_2d_block(dx=tank_spacing,
                                length=(tank_layers - 1) * tank_spacing,
                                height=tank_height,
                                center=[0., 0.])
    xleft += min(xf) - max(xleft) - tank_spacing
    yleft += min(yf) - min(yleft)

    xright = xleft + abs(min(xleft)) + tank_length + tank_spacing
    yright = yleft

    xbottom, ybottom = get_2d_block(dx=tank_spacing,
                                    length=max(xright) - min(xleft),
                                    height=(tank_layers - 1) * tank_spacing,
                                    center=[0., 0.])
    xbottom += min(xleft) - min(xbottom)
    ybottom += min(yleft) - max(ybottom) - tank_spacing

    xtop = np.array([])
    ytop = np.array([])
    if close is True:
        xtop, ytop = get_2d_block(dx=tank_spacing,
                                  length=max(xright) - min(xleft),
                                  height=(tank_layers - 1) * tank_spacing,
                                  center=[0., 0.])
        xtop += min(xleft) - min(xtop)
        ytop += max(yleft) - min(ytop) - tank_spacing

    x = np.concatenate([xleft, xright, xbottom, xtop])
    y = np.concatenate([yleft, yright, ybottom, ytop])

    return x, y
