import overlapadd2 as ola
import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline
radii = np.fromfile('r.bin', np.float64)
hyper = np.fromfile('hyper.bin', np.float64)
interpHyper = InterpolatedUnivariateSpline(radii, hyper, ext='raise')

import halfbandfilter as hb
from scipy.signal import convolve2d

import gdal, gdalconst
import textureShading as tex

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

def memmapInput(binfilename, memmapMode='r', dtype=np.int16, shape=None):
    if shape is None:
        fileHandle = gdal.Open(binfilename, gdalconst.GA_ReadOnly)
        width, height = (fileHandle.RasterXSize, fileHandle.RasterYSize)
        shape = (height, width)

    inmap = np.memmap(binfilename, dtype=dtype, mode=memmapMode,
            shape=shape)
    return inmap
    
def memmapOutput(fname, shape):
    return np.memmap(fname, dtype=np.float32, mode='w+', shape=tuple(shape))

def run(elevBinName, hankelTaps=960, L=(3500, 3500), verbose=True, workingDir='./'):
    elevation = memmapInput(elevBinName)

    hankelFilter = makeHalfHankel(hankelTaps)
    if verbose:
        print 'filter size:', hankelFilter.shape
    outSize = np.array(elevation.shape) + hankelFilter.shape - 1

    texture = memmapOutput(workingDir + 'tex.bin', outSize)
    texture = ola.overlapadd2(elevation, hankelFilter, y=texture, L=L, verbose=True)
    if verbose:
        print("Done with filtering")

    origSize = elevation.shape
    filterSize = hankelFilter.shape
    subtexture = memmapOutput(workingDir + 'subtex.bin', elevation.shape)
    subtexture[:] = texture[filterSize[0]/2 : filterSize[0]/2+origSize[0], 
                            filterSize[1]/2 : filterSize[1]/2+origSize[1] ][:]
    subtexture.flush()
    if verbose:
        print("Done with extracting sub-texture")

    return subtexture

class RunningHist:
    def __init__(self, lo, hi, nbins=1000):
        self.range = (lo, hi)
        self.nbins = nbins
        self.bins = np.linspace(lo, hi, nbins+1)
        self.frequencies = np.zeros(nbins)
    def __call__(self, data):
        newFreq, _ = np.histogram(data, bins=self.bins, range=self.range, density=False)
        self.frequencies += newFreq
    def pmf(self):
        return (self.frequencies / (np.sum(self.frequencies) * np.diff(self.bins[:2])), self.bins)

def bigMinmax(x, mask):
    lo = np.finfo(x.dtype).max
    hi = np.finfo(x.dtype).min
    for i in xrange(len(x)):
        row = x[i, mask[i].astype(np.bool)]
        if row.size > 0:
            lo = min(lo, row.min())
            hi = max(hi, row.max())
    return [lo, hi]

def estimateLimits(subtexture, mask, percentiles, postPercentileScale, 
                   skip=10, nbins=1000, small=False):
    if small:
        hist, binEdges = np.histogram(subtexture[mask][::skip],
                bins=nbins, density=True)
    else:
        lohi = bigMinmax(subtexture, mask)
        h = RunningHist(*lohi, nbins=1000)
        for i in xrange(len(subtexture)):
            h(subtexture[i, mask[i].astype(np.bool)][::skip])
        hist, binEdges = h.pmf()
    distribution = np.cumsum(hist) * np.diff(binEdges[:2])
    
    loIdx = np.sum(distribution < percentiles[0] / 100.0)
    hiIdx = np.sum(distribution < percentiles[1] / 100.0)
    
    if loIdx == 0 or hiIdx == nbins or hiIdx == loIdx:
        print("Histogram with {} bins insufficient to estimate percentiles: lo={} and hi={}".format(nbins, loIdx, hiIdx))

    rawLimits = np.array([binEdges[loIdx], binEdges[hiIdx]])
    limits = rawLimits * postPercentileScale
    return tuple(limits)

def postProcess(subtexture, maskName, limits=None, 
                percentiles=[0.1, 99.9], postPercentileScale=[1.0, 1.0],
                texTifName=None, verbose=True):
    # Find no-data mask from original file, which should be the water portions
    # of the world according to Natural Earth vector data
    landMask = memmapInput(maskName, memmapMode='r+', dtype=np.uint8)
    if verbose:
        print("Generated no-data mask")
    
    if limits is None:
        # Get a soft-estimate of the color limits within `percentiles`
        # percentiles
        limits = estimateLimits(subtexture, landMask, percentiles,
                postPercentileScale)
        if verbose:
            print("Estimated fuzzy limits", limits)

    # Within those limits, rescale to 1-255 (leaving 0 as a new no-data value)
    landMask = tex.touintChunked(subtexture, *limits, dtype=np.uint8, ret=landMask)
    if verbose:
        print("Converted to byte image")

    if texTifName is None:
        return limits

    # In the byte image, the no-data value is 0. Don't reinstate it here: do so outside Python
    newNDV = 0

    # Save!
    driver = gdal.GetDriverByName('GTiff')
    NDV, xsize, ysize, GeoT, Projection, DataType = tex.GetGeoInfo(maskName)

    tex.CreateGeoTiff(texTifName, landMask, driver, newNDV, xsize, ysize, 
            GeoT, Projection, gdalconst.GDT_Byte)
    if verbose:
        print("Done saving.")
    return limits

def disp(tex):
    import pylab
    pylab.ion()
    pylab.figure()
    pylab.imshow(tex)
    pylab.colorbar()
    pylab.set_cmap('gray')
    quantiles = [0.5, 99.9];
    pylab.clim(np.percentile(t, quantiles))

"""
Starting with a GeoTIF file, run:
$ gdal_calc.py -A INPUT.tif --outfile=OUTPUT.bin --calc="A*(A>-30000)" --NoDataValue=0 --format=ENVI
This converts the GeoTIF file to a flat binary file (the ENVI standard), and
replaces the original NoDataValue of -32K with 0, so we don't have to do it
here.

Then call 
```
run("OUTPUT.bin", hankelTaps=1984, L=(7000, 7000))
```
The filter will be of size `hankelTaps / 2 + 32`, or, in the above example, 1024
by 1024. `L` is a two-tuple and gives the size of the sub-convolutions that
overlap-add will perform using the FFT. Make this as big as possible, given your
memory limits. It doesn't need to be square. You want to set `L` such that `L +
(size of filter) - 1` is just below a power-of-two. In this example, `7000 +
1024 - 1` is 8023, which is just below 8K, or 8192.
"""
if __name__ == '__main__':
    datadir = '/Users/ahmed.fasih/Documents/Personal/textureShading-Asia/Data/'

    if not True:  # Testing
        print("Test")
        intif = datadir + 'east_0.1.tif'
        inbin = datadir + 'east_0.1.bin'
        t = run(inbin, hankelTaps=288, L=(500, 500))
    elif False: # 250 meter data
        print("250 meter data, east or west")
        intif = datadir + 'east-land.tif'
        inbin = datadir + 'east-land.bin'

        intif = datadir + 'w-land.tif'
        inbin = datadir + 'w-land.bin'
        t = run(inbin, hankelTaps=4032, L=(6100, 14300))
    else:  # 90 meter data!
        print("90 meter data!")
        datadir = '/Users/ahmed.fasih/Documents/Personal/textureShading-Asia/'
        wdir = '/Volumes/SeagateBack/Fasih/90m/'
        inbin = datadir + 'land.bin'
        run(inbin, hankelTaps=12224, L=(2000, 10000), workingDir=wdir)
    
    # postProcess(t, intif)

def makeTileCommand(filename):
    ll = tex.filenameToLatsLons(filename)
    NDV, xsize, ysize, GeoT, Projection, DataType = tex.GetGeoInfo(filename)
    return "gdal_translate -of VRT -a_srs EPSG:4326 -gcp 0 0 {} {} -gcp {} 0 {} {} -gcp {} {} {} {} {} ne.vrt && gdalwarp -of VRT -t_srs EPSG:4326 ne.vrt ne2.vrt && python gdal2tiles_parallel.py -p mercator -a 0,0,0 --config GDAL_CACHEMAX 1024 ne2.vrt".format(ll['lonExt'][0], ll['latExt'][1], xsize, ll['lonExt'][1], ll['latExt'][1], xsize, ysize, ll['lonExt'][1], ll['latExt'][0], filename)
    

"""
To remove lakes *in-place*:

$ gdal_rasterize -sql "SELECT * FROM ne_10m_lakes WHERE scalerank=0" -burn 0  Data/ne_10m_lakes/ne_10m_lakes.shp INPUT.tif
"""

def sph2cart(az, el, r=6371e3):
    x = r * np.cos(el) * np.cos(az)
    y = r * np.cos(el) * np.sin(az)
    z = r * np.sin(el)
    return (x, y, z)

def sphd2cart(az, el, r=6371e3):
    return sph2cart(np.deg2rad(az), np.deg2rad(el), r)

"""
import numpy as np
import gdal, gdalconst
import textureShading as tex
import big

mask = big.memmapInput('land-mask.bin',dtype=np.uint8)
text = big.memmapInput('/Volumes/SeagateBack/Fasih/90m/subtex.bin', shape=mask.shape, dtype=np.single)
big.postProcess(text, 'land-mask.bin')


mask = big.memmapInput('land-mask-0.1.bin',dtype=np.uint8)
text = big.memmapInput('subtex.bin',shape=mask.shape,dtype=np.single)
big.postProcess(text, 'land-mask-0.1.bin', texTifName='small-runninghist.tif')
"""
