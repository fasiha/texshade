import numpy as np
from scipy.signal import decimate
import gdal, gdalconst


def touint16(x, cmin, cmax):
    from scipy.interpolate import interp1d
    interpolator = interp1d([cmin, cmax], [0, 2 ** 16 - 1 - 1e-6])
    clamp = lambda x, lo, hi: np.minimum(np.maximum(x, lo), hi)
    y = interpolator(clamp(x, cmin, cmax)).astype(np.uint16)
    return y


texture = np.load('texpy.npy')
if True:
    filename = 'Data/east_0.1.tif'
    fileHandle = gdal.Open(filename, gdalconst.GA_ReadOnly)
    band = fileHandle.GetRasterBand(1)
    orig = band.ReadAsArray()
    texture[orig < 0] = 0
    orig = [] # clear

colorlimits = np.percentile(texture.ravel(), [1.0 / 2, 99.9])
newlimits = colorlimits * [.55, .4]
newlimits = colorlimits * [.25, .2]

save = True
if save:
    import png
    png.from_array(touint16(texture, *colorlimits), 'L').save('pytex16.png')

show = False
if show:
    import pylab
    pylab.ion()

    stage = lambda x, ax: decimate(x, 2, n=2, axis=ax)
    horiz = lambda x: stage(stage(stage(x, 1), 1), 1)
    vert = lambda x: stage(stage(stage(x, 0), 0), 0)
    down = horiz(vert(texture))

    pylab.imshow(down, cmap=pylab.cm.gray)
    pylab.colorbar()

    colorlimits = np.percentile(texture.ravel(), [1.0 / 2, 99.9])
    pylab.clim(newlimits)
