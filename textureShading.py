import numpy as np
import gdal, gdalconst
import numpy.fft as fft
import scipy.fftpack as scifft
from osgeo import osr

nextpow2 = lambda v: 2 ** int(np.ceil(np.log2(v)))


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
    H2 = None # potentially trigger GC here to reclaim H2's memory
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
    latSpacing = np.diff(latExt)[0] / fileHandle.RasterYSize

    lons = np.linspace(lonExt.min(), lonExt.max(), fileHandle.RasterXSize + 1)
    lonSpacing = np.diff(lonExt)[0] / fileHandle.RasterXSize
    return {'lats' : lats, 'lons' : lons, 'latSpacing' : latSpacing,
            'lonSpacing' : lonSpacing, 'latExt' : latExt, 'lonExt' : lonExt}

def scaleByPixelArea(x, lats, lons):
    areas = latLonRectToArea(lats, lons)
    return x / (areas[:, np.newaxis] / areas.max()).astype(x.dtype)

def filenameToNoDataMask(fname, bandId = 1):
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)
    band = fileHandle.GetRasterBand(bandId)
    return band.ReadAsArray() == band.GetNoDataValue()

def filenameToNoDataValue(fname, bandId = 1):
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)
    band = fileHandle.GetRasterBand(bandId)
    return band.GetNoDataValue()

def filenameToData(fname, dtype = np.float32):
    """Reads all bands"""
    fileHandle = gdal.Open(fname, gdalconst.GA_ReadOnly)
    result = np.squeeze(np.dstack([fileHandle.GetRasterBand(n+1).ReadAsArray() for
        n in range(fileHandle.RasterCount)]))
    if dtype is not None:
        return result.astype(dtype)
    return result

def filenameToTexture(fname, alpha = 0.5, ndvReplacement = 0.0, verbose = True):
    ndv = filenameToNoDataValue(fname) # no data value
    latsLons = filenameToLatsLons(fname)
    x = filenameToData(fname)
    if verbose: print "Loaded data"

    x[x==ndv] = ndvReplacement
    x = scaleByPixelArea(x, latsLons['lats'], [0, latsLons['lonSpacing']])
    if verbose: print "Done preprocessing"
    
    x = sciMoreReal(x, alpha)
    if verbose: print "Computed texture"
    
    return x

def touintChunked(x, cmin, cmax, dtype=np.uint16, verbose=False, ret=None,
        mask=None):
    maxval = dtype(2 ** (8 * dtype().itemsize) - 1e-3)
    slope = (maxval - 1.0) / (cmax - cmin)

    # allocate up-front if needed
    if ret is None:
        ret = np.zeros(x.shape, dtype=dtype)
    
    # clamp x between cmin and cmax. Map [cmin, cmax] to [1, 2**depth-1-1e-3]
    # linearly. Reserve 0 for nodata.
    for i in xrange(x.shape[0]):
        if verbose: print("Row {} of {}".format(i, x.shape[0]))
        if mask is None:
            thisMask = np.newaxis
        else:
            thisMask = mask[i].astype(np.bool)
        ret[i, thisMask] = ((x[i, thisMask] <= cmin) * 1 + 
                (x[i, thisMask] >= cmax) * maxval + 
                np.logical_and(x[i, thisMask] > cmin, 
                               x[i, thisMask] < cmax) * 
                    (slope * (x[i] - cmin) + 1))
    return ret

def touint(x, cmin, cmax, dtype = np.uint16):
    # clamp x between cmin and cmax
    x[x < cmin] = cmin
    x[x > cmax] = cmax
    # map [cmin, cmax] to [1, 2**depth-1-1e-6] linearly
    maxval = 2 ** (8 * dtype().itemsize) - 1e-6
    slope = (maxval - 1.0) / (cmax - cmin)
    ret = slope * (x - cmin) + 1
    return (ret).astype(dtype)

def dataToPercentileLimits(texture, percentiles, postPercentileScale):
    colorlimits = np.percentile(texture.ravel(), percentiles)
    return colorlimits * postPercentileScale

def dataToPng16(data, fname = 'pytex16.png', percentiles = [0.5, 99.9],
        postPercentileScale = [1.0, 1.0]):
    import png
    limits = dataToPercentileLimits(data, percentiles, postPercentileScale)
    png.from_array(touint(data, *limits, depth = 16), 'L').save(fname)

def geoFileToStruct(fname):
    proj = osr.SpatialReference()

    source = gdal.Open(fname, gdalconst.GA_ReadOnly)
    proj.ImportFromWkt(source.GetProjectionRef())

    band = source.GetRasterBand(1)
    
    ret = dict(ndv=band.GetNoDataValue(),
            width=source.RasterXSize,
            height=source.RasterYSize,
            transform=source.GetGeoTransform(),
            projection=proj,
            dtype=band.DataType)
    return ret

def structToGeoFile(fname, params, array=None, driverString='GTiff'):
    driver = gdal.GetDriverByName(driverString)
    data = driver.Create(fname, params['width'], params['height'], 1, params['dtype'])
    data.SetGeoTransform(params['transform'])
    data.SetProjection(params['projection'].ExportToWkt())
    if params['ndv'] is not None:
        data.GetRasterBand(1).SetNoDataValue(params['ndv'])
    if array is not None:
        data.GetRasterBand(1).WriteArray(array)
    return data

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
                  xsize, ysize, GeoT, Projection, DataType, numBands=1):
    # Set up the dataset
    DataSet = driver.Create( Name, xsize, ysize, numBands, DataType, options = [ 'COMPRESS=LZW' ] )
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection( Projection.ExportToWkt() )
    # Write the array
    if numBands == 1:
        DataSet.GetRasterBand(numBands).WriteArray( Array )
        if NDV is not None: DataSet.GetRasterBand(numBands).SetNoDataValue(NDV)
    else:
        for bid in range(numBands):
            DataSet.GetRasterBand(bid + 1).WriteArray( Array[:,:,bid] )
            if NDV is not None: DataSet.GetRasterBand(bid + 1).SetNoDataValue(NDV)
    return Name

def dataToGTiff(data, infile, outfile='pytex.tiff', depth=16, 
        percentiles=[0.5, 99.9], postPercentileScale=[1.0, 1.0], sample=1.0,
        remask=True):
    if depth == 16:
        gdalDataType = gdalconst.GDT_UInt16
        npDataType = np.uint16
    elif depth == 32:
        gdalDataType = gdalconst.GDT_UInt32
        npDataType = np.uint32
    elif depth == 8:
        gdalDataType = gdalconst.GDT_Byte
        npDataType = np.uint8
    else:
        raise ValueError('Only 8-, 16-, or 32-bit GTiffs allowed')
    sample = max(min(sample, 1.0), 1e-6)
    skip = int(1.0 / sample)
    if data.size <= skip:
        raise ValueError('Sample rate too low, less than two samples will be checked')

    driver = gdal.GetDriverByName('GTiff')
    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(infile)

    # Calculate percentiles of non-masked values
    if remask:
        ndvMask = filenameToNoDataMask(infile)
        limits = dataToPercentileLimits(data[np.logical_not(ndvMask)][::skip],
                percentiles, postPercentileScale)
    else:
        limits = dataToPercentileLimits(data[::skip, ::skip], percentiles,
                postPercentileScale)

    # Remap to output format
    data = touint(data, *limits, dtype = npDataType)
    # print "Chose limits", limits, "minmax of result", [np.min(data), np.max(data)]

    # Reinstate mask
    newNDV = 0.0
    data[ndvMask] = newNDV
    
    CreateGeoTiff(outfile, data, driver, newNDV, xsize, ysize, GeoT, Projection,
            gdalDataType)

