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


filename = 'Data/east_0.1.tif'
fileHandle = gdal.Open(filename, gdalconst.GA_ReadOnly)

ext = np.array(GetExtent(fileHandle.GetGeoTransform(), fileHandle.RasterXSize,
        fileHandle.RasterYSize))
latExt = np.unique(ext[:, 1])
lonExt = np.unique(ext[:, 0])
lats = np.linspace(latExt.max(), latExt.min(), fileHandle.RasterYSize + 1)
lons = [0, np.diff(lonExt) / fileHandle.RasterXSize]
areas = latLonRectToArea(lats, lons)

band = fileHandle.GetRasterBand(1)
data = band.ReadAsArray().astype(np.float32)
if not False:
    # land only
    data = np.maximum(0, data)
    # scale by area: polar areas have "higher" elevation than equatorial areas
    data = data / (areas[:,np.newaxis] / areas.max())
print "Read data"

data = sciMoreReal(data, 0.5)
print "Computed texture"

np.save('texpy.npy', data)
print "Saved and exiting"
