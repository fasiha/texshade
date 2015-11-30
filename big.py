import overlapadd2 as ola
import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline
radii = np.fromfile('r.bin', np.float64)
hyper = np.fromfile('hyper.bin', np.float64)
interpHyper = InterpolatedUnivariateSpline(radii, hyper, ext='raise')

import halfbandfilter as hb
from scipy.signal import convolve2d
import gdal, gdalconst

def vec(v):
    return v.reshape(v.size, -1)

def makeAxis(n):
    return np.arange(-n/2.0, n/2.0)

def makeRadii(N):
    return np.sqrt(vec(makeAxis(N)) ** 2 + vec(makeAxis(N)).T ** 2)

def makeHankel(N):
    r = makeRadii(N)
    return interpHyper(r)

def makeHalfHankel(N, hbTaps=64):

    spatialFilter = makeHankel(N)
    hbFilter = hb.design(hbTaps)
    doubleFilter = convolve2d(convolve2d(spatialFilter, vec(hbFilter)),
            vec(hbFilter).T)

    if N % 4 == 0:
        finalFilter = doubleFilter[:-1:2, :-1:2]
    else:
        finalFilter = doubleFilter[1:-1:2, 1:-1:2]
    
    return finalFilter

def memmapInput(origfilename, binfilename, memmapMode='c'):
    fileHandle = gdal.Open(origfilename, gdalconst.GA_ReadOnly)
    width, height = (fileHandle.RasterXSize, fileHandle.RasterYSize)

    inmap = np.memmap(binfilename, dtype=np.float32, mode=memmapMode,
            shape=(height, width))
    return inmap
    
def memmapOutput(fname, shape):
    return np.memmap(fname, dtype=np.float32, mode='w+', shape=tuple(shape))

def run(origfilename, binfilename, hankelTaps=960, L=(3500, 3500), verbose=True):
    elevation = memmapInput(origfilename, binfilename)
    elevation[elevation < 0] = 0.0

    hankelFilter = makeHalfHankel(hankelTaps)
    if verbose:
        print 'filter size:', hankelFilter.shape
    outSize = np.array(elevation.shape) + hankelFilter.shape - 1

    texture = memmapOutput('tex.bin', outSize)

    texture = ola.overlapadd2(elevation, hankelFilter, y=texture, L=L, verbose=True)

    return texture

def disp(tex):
    import pylab
    pylab.ion()
    pylab.figure()
    pylab.imshow(tex)
    pylab.colorbar()
    pylab.set_cmap('gray')
    quantiles = [0.5, 99.9];
    pylab.clim(np.percentile(t, quantiles))


if __name__ == '__main__':
# t = run('chgis_dem.tif', 'foo.bin', 960)
    t = run('/Users/ahmed.fasih/Documents/Personal/textureShading-Asia/Data/SRTM_NE_250m.tif',
        '/Users/ahmed.fasih/Documents/Personal/textureShading-Asia/Data/ne-land-single.bin',
        960, L=(2*3500, 2*3500))

