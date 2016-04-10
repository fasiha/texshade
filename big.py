import gdal, gdalconst
import textureShading as tex

import itertools as itools
import overlapadd2 as ola
import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline
radii = np.fromfile('r.bin', np.float64)
hyper = np.fromfile('hyper.bin', np.float64)
interpHyper = InterpolatedUnivariateSpline(radii, hyper, ext='raise')

import halfbandfilter as hb
from scipy.signal import convolve2d

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

def mkchunks(total, chunkSize):
    numChunks = np.ceil(float(total) / chunkSize).astype(long)
    idxToRange = lambda i: xrange(i * chunkSize, min((i + 1) * chunkSize, total))
    return itools.imap(idxToRange, xrange(numChunks))

def mkchunksN(totals, chunkSizes):
    if len(chunkSizes) == 1:
        chunkSizes *= len(totals)

    iterators = itools.imap(mkchunks, totals, chunkSizes)
    return itools.product(*iterators)

def geoFileToChunk(fname, start0, end0, start1, end1,bandnum=1):
    src = gdal.Open(fname, gdalconst.GA_ReadOnly)
    band = src.GetRasterBand(bandnum)
    return band.ReadAsArray(start1, start0, end1-start1, end0-start0)

def run(elevBinName, hankelTaps=960, L=(3500, 3500), verbose=True,
        workingDir='./'):
    # Hankel filter
    hankelFilter = makeHalfHankel(hankelTaps)
    if verbose:
        print 'filter size:', hankelFilter.shape
    filterSize = hankelFilter.shape

    # Prep input data
    inHandle = gdal.Open(elevBinName, gdalconst.GA_ReadOnly)
    inBand = inHandle.GetRasterBand(1)
    bandNDV = inBand.GetNoDataValue()

    origSize = (inHandle.RasterYSize, inHandle.RasterXSize)
    outSize = np.array(origSize) + filterSize - 1

    cleaner = lambda arr: arr * np.logical_not(np.isclose(arr, bandNDV))
    getElevation = lambda *args: geoFileToChunk(elevBinName, *args)
    getCleanElevation = lambda *args: cleaner(getElevation(*args))

    # Prep output
    outputName = workingDir + 'tex-2.bin'
    outParams = tex.geoFileToStruct(elevBinName)
    outParams['dtype'] = gdal.GetDataTypeByName('float32')
    outParams['width'] = outSize[1]
    outParams['height'] = outSize[0]
    # For details about 'transform' list here, see:
    #   http://www.gdal.org/classGDALDataset.html#af9593cc241e7d140f5f3c4798a43a668
    affine = list(outParams['transform'])
    affine[0] += affine[1] * -filterSize[1]/2
    affine[3] += affine[5] * -filterSize[0]/2
    outParams['transform'] = tuple(affine)
    createHandle = tex.structToGeoFile(outputName, outParams, driverString='ENVI')
    createHandle.FlushCache()
    del createHandle  # otherwise, GDAL (and maybe memmap?) can't open file below


    # Produce output
    ndv = outParams['ndv']

    texture = memmapOutput(outputName, outSize)

    def textureAdderClean(start0,end0,start1,end1,arr):
        inr0 = xrange(start0, end0)
        inr1 = xrange(start1, end1)

        outr0 = np.array(inr0) - filterSize[0] / 2
        outr1 = np.array(inr1) - filterSize[1] / 2

        good0 = np.logical_and(outr0 >= 0, outr0 < origSize[0])
        good1 = np.logical_and(outr1 >= 0, outr1 < origSize[1])

        el0_start = outr0[good0][0]
        el0_end = outr0[good0][-1]+1
        el1_start = outr1[good1][0]
        el1_end = outr1[good1][-1]+1

        ndvmask = np.ones(arr.shape, dtype=bool)
        ndvmask[np.c_[np.where(good0)], good1] = np.logical_not(np.isclose(
            ndv, getElevation(el0_start, el0_end, el1_start, el1_end)))

        texture[start0:end0, start1:end1] = (np.logical_not(ndvmask) * ndv +
                ndvmask * (arr + texture[start0:end0, start1:end1]))
    def textureAdder(start0,end0,start1,end1,arr):
        texture[start0:end0, start1:end1] += arr

    ola.overlapadd2(getCleanElevation, hankelFilter, y=textureAdderClean, L=L,
            verbose=True, Na=origSize)
    texture.flush()
    if verbose:
        print("Done with filtering")


    # Clean up texture output: burn NDVs into FIR-induced border and whever
    # original had NDV
    outHandle = gdal.Open(outputName, gdalconst.GA_Update)
    outBand = outHandle.GetRasterBand(1)
    vertSlab = ndv * np.ones([outSize[0], filterSize[1]/2], dtype=float)
    horizSlab = ndv * np.ones([filterSize[1]/2, outSize[1]], dtype=float)
    outBand.WriteArray(horizSlab)
    outBand.WriteArray(vertSlab)
    # right
    outBand.WriteArray(vertSlab[:,:-1], xoff=filterSize[1]/2 + origSize[1], yoff=0)
    # bottom
    outBand.WriteArray(horizSlab[:-1,:], xoff=0, yoff=filterSize[0]/2 + origSize[0])

    return texture

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
    import sys
    if len(sys.argv) == 1:
        setup = "test"
    else:
        setup = sys.argv[-1]

    datadir = '/Users/ahmed.fasih/Downloads/gdal-demo-cgiar-srtm/'

    if setup == 'test':  # Testing
        nedtif = '/Users/ahmed.fasih/Documents/Personal/texshading/ned-small.tif'
        srtmtif = '/Users/ahmed.fasih/Downloads/gdal-demo-cgiar-srtm/SRTM-small.tif'
        intif = nedtif or srtmtif
        t = run(intif, hankelTaps=288, L=(500, 500))
    elif setup == 'cgiar-250': # 250 meter data
        print("250 meter data, east or west")
        intif = datadir + 'east-land.tif'
        inbin = datadir + 'east-land.bin'

        intif = datadir + 'w-land.tif'
        inbin = datadir + 'w-land.bin'
        t = run(inbin, hankelTaps=4032, L=(6100, 14300))
    elif setup == 'cgiar':  # 90 meter data!
        print("90 meter data!")
        datadir = '/Users/ahmed.fasih/Documents/Personal/textureShading-Asia/'
        wdir = '/Volumes/SeagateBack/Fasih/90m/'
        inbin = datadir + 'land.bin'
        run(inbin, hankelTaps=12224, L=(2000, 10000), workingDir=wdir)
    elif setup == 'ned': # 10m NED data!
        t = run('/srv/data/fasih/NED/vrt/ned.vrt', hankelTaps=12224, L=(10000,
            25000), workingDir='/srv/data/fasih/NED/vrt/')

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
