import numpy as np
import gdal, gdalconst
import numpy.fft as fft
import scipy.fftpack as scifft
from osgeo import osr

nextpow2 = lambda v: 2 ** np.ceil(np.log2(v)).astype(int)


def sciMoreReal(x, alpha, verbose = True):
    Nyx = nextpow2(x.shape)

    fy = scifft.rfftfreq(Nyx[0])[:, np.newaxis].astype(x.dtype)
    fx = scifft.rfftfreq(Nyx[1])[np.newaxis, :].astype(x.dtype)
    H2 = np.sqrt(fx ** 2 + fy ** 2) ** alpha
    if verbose: print "Generated filter"

    rfft2 = lambda x: scifft.rfft(scifft.rfft(x, Nyx[1], 1, True), Nyx[0], 0, True)
    irfft2 = lambda X: scifft.irfft(scifft.irfft(X, axis=0, overwrite_x=True), overwrite_x=True)

    xr = rfft2(x) * H2
    if verbose: print "Completed frequency domain operations"
    H2 = [] # potentially trigger GC here to reclaim H2's memory
    xr = irfft2(xr)
    if verbose: print "Back to spatial-domain"

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

def filenameToNoDataMask(fname, bandId = 1):
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)
    band = fileHandle.GetRasterBand(bandId)
    return band.ReadAsArray() == band.GetNoDataValue()

def filenameToNoDataValue(fname, bandId = 1):
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)
    band = fileHandle.GetRasterBand(bandId)
    return band.GetNoDataValue()

def filenameToData(fname, dtype = np.float32, bandId = 1):
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)
    band = fileHandle.GetRasterBand(bandId)
    return band.ReadAsArray().astype(dtype)

def filenameToTexture(fname, alpha = 0.5, ndvReplacement = 0.0, verbose = True):
    ndv = filenameToNoDataValue(fname) # no data value
    latsLons = filenameToLatsLons(fname)
    x = filenameToData(fname)
    if verbose: print "Loaded data"

    x[x==ndv] = ndvReplacement
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

def dataToPng16(data, fname = 'pytex16.png', percentiles = [0.5, 99.9],
        postPercentileScale = [1.0, 1.0]):
    import png
    limits = dataToPercentileLimits(data, percentiles, postPercentileScale)
    png.from_array(touint16(data, *limits), 'L').save(fname)

# Adapted from EddyTheB, http://gis.stackexchange.com/a/57006/8623
def GetGeoInfo(FileName):
    SourceDS = gdal.Open(FileName, gdalconst.GA_ReadOnly)
    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())
    DataType = SourceDS.GetRasterBand(1).DataType
    DataType = gdal.GetDataTypeName(DataType)
    return NDV, xsize, ysize, GeoT, Projection, DataType

def CreateGeoTiff(Name, Array, driver, NDV, 
                  xsize, ysize, GeoT, Projection, DataType):
    if DataType == 'Float32':
        DataType = gdal.GDT_Float32
    NewFileName = Name+'.tif'
    # Set nans to the original No Data Value
    if NDV is not None: Array[np.isnan(Array)] = NDV
    # Set up the dataset
    DataSet = driver.Create( NewFileName, xsize, ysize, 1, DataType )
            # the '1' is for band 1.
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection( Projection.ExportToWkt() )
    # Write the array
    DataSet.GetRasterBand(1).WriteArray( Array )
    if NDV is not None: DataSet.GetRasterBand(1).SetNoDataValue(NDV)
    return NewFileName

def dataToGTiff16(data, infile, outfile = 'pytex16.tiff', percentiles = [0.5, 99.9],
        postPercentileScale = [1.0, 1.0]):
    driver = gdal.GetDriverByName('GTiff')
    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(infile)

    limits = dataToPercentileLimits(data, percentiles, postPercentileScale)
    data = touint16(data, *limits)
    
    CreateGeoTiff(outfile, data, driver, NDV, xsize, ysize, GeoT, Projection,
            gdalconst.GDT_UInt16)

