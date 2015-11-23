import numpy as np
import gdal, gdalconst
import numpy.fft as fft
import scipy.fftpack as scifft

nextpow2 = lambda v: 2 ** np.ceil(np.log2(v)).astype(int)


def sciMoreReal(x, alpha):
    Nyx = nextpow2(x.shape)

    fy = scifft.rfftfreq(Nyx[0])[:, np.newaxis].astype(x.dtype)
    fx = scifft.rfftfreq(Nyx[1])[np.newaxis, :].astype(x.dtype)
    H2 = np.sqrt(fx ** 2 + fy ** 2) ** alpha

    rfft2 = lambda x: scifft.rfft(scifft.rfft(x, Nyx[1], 1, True), Nyx[0], 0, True)
    irfft2 = lambda X: scifft.irfft(scifft.irfft(X, axis=0, overwrite_x=True), overwrite_x=True)

    xr = rfft2(x) * H2
    H2 = [] # potentially trigger GC here to reclaim H2's memory
    xr = irfft2(xr)

    return xr[:x.shape[0], :x.shape[1]]

# From http://gis.stackexchange.com/a/57837/8623
def GetExtent(gt,cols,rows):
    ''' Return list of corner coordinates from a geotransform

        @type gt:   C{tuple/list}
        @param gt: geotransform
        @type cols:   C{int}
        @param cols: number of columns in the dataset
        @type rows:   C{int}
        @param rows: number of rows in the dataset
        @rtype:    C{[float,...,float]}
        @return:   coordinates of each corner
    '''
    ext=[]
    xarr=[0,cols]
    yarr=[0,rows]

    for px in xarr:
        for py in yarr:
            x=gt[0]+(px*gt[1])+(py*gt[2])
            y=gt[3]+(px*gt[4])+(py*gt[5])
            ext.append([x,y])
        yarr.reverse()
    return ext

def ReprojectCoords(coords,src_srs,tgt_srs):
    ''' Reproject a list of x,y coordinates.

        @type geom:     C{tuple/list}
        @param geom:    List of [[x,y],...[x,y]] coordinates
        @type src_srs:  C{osr.SpatialReference}
        @param src_srs: OSR SpatialReference object
        @type tgt_srs:  C{osr.SpatialReference}
        @param tgt_srs: OSR SpatialReference object
        @rtype:         C{tuple/list}
        @return:        List of transformed [[x,y],...[x,y]] coordinates
    '''
    trans_coords=[]
    transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords

def latLonRectToArea(lats, lons, radius = 6371.0):
    """lats, lons : array-like in degrees"""
    lats = np.array(lats) * np.pi / 180.0
    return np.pi / 180.0 * radius ** 2 * np.abs(np.diff(np.sin(lats)) * np.diff(lons))

def landOnly(x, cutoff = 0.0):
    return np.maximum(cutoff, x)

def filenameToLatsLons(fname):
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)

    ext = np.array(GetExtent(fileHandle.GetGeoTransform(), fileHandle.RasterXSize,
            fileHandle.RasterYSize))
    latExt = np.unique(ext[:, 1])
    lonExt = np.unique(ext[:, 0])
    lats = np.linspace(latExt.max(), latExt.min(), fileHandle.RasterYSize + 1)
    lons = [0, np.diff(lonExt) / fileHandle.RasterXSize]
    return {'lats' : lats, 'lons' : lons}

def scaleByPixelArea(x, lats, lons):
    areas = latLonRectToArea(lats, lons)
    return x / (areas[:, np.newaxis] / areas.max())

def filenameToLandMask(fname, cutoff = 0, bandId = 1):
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)
    band = fileHandle.GetRasterBand(bandId)
    return band.ReadAsArray() < cutoff

def filenameToData(fname, dtype = np.float32, bandId = 1):
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)
    band = fileHandle.GetRasterBand(bandId)
    return band.ReadAsArray().astype(dtype)

def filenameToTexture(fname, alpha = 0.5, clipNegative = True, verbose = True):
    latsLons = filenameToLatsLons(fname)
    x = filenameToData(fname)
    if verbose: print "Loaded data"
    if clipNegative:
        x = landOnly(x)
    x = scaleByPixelArea(x, latsLons['lats'], latsLons['lons'])
    if verbose: print "Done preprocessing"
    x = sciMoreReal(x, alpha)
    if verbose: print "Computed texture"
    return x

def touint16(x, cmin, cmax):
    from scipy.interpolate import interp1d
    interpolator = interp1d([cmin, cmax], [0, 2 ** 16 - 1 - 1e-6])
    clamp = lambda x, lo, hi: np.minimum(np.maximum(x, lo), hi)
    y = interpolator(clamp(x, cmin, cmax)).astype(np.uint16)
    return y

def dataToPercentileLimits(texture, percentiles, postPercentileScale):
    colorlimits = np.percentile(texture.ravel(), percentiles)
    return colorlimits * postPercentileScale

def dataToPng(data, fname = 'pytex16.png', percentiles = [0.5, 99.9],
        postPercentileScale = [1.0, 1.0]):
    import png
    limits = dataToPercentileLimits(data, percentiles, postPercentileScale)
    png.from_array(touint16(data, *limits), 'L').save(fname)


