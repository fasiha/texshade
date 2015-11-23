# Adapted from EddyTheB, http://gis.stackexchange.com/a/57006/8623
import numpy as np
import gdal
from gdalconst import *
from osgeo import osr

# Function to read the original file's projection:
def GetGeoInfo(FileName):
    SourceDS = gdal.Open(FileName, GA_ReadOnly)
    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())
    DataType = SourceDS.GetRasterBand(1).DataType
    DataType = gdal.GetDataTypeName(DataType)
    return NDV, xsize, ysize, GeoT, Projection, DataType

# Function to write a new file.
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

# Open the original file
FileName = '/Users/fasih/Dropbox/Data/cgiar/0.01_small_East.tif'
DataSet = gdal.Open(FileName, GA_ReadOnly)
# Get the first (and only) band.
Band = DataSet.GetRasterBand(1)
# Open as an array.
Array = Band.ReadAsArray()
# Get the No Data Value
NDV = Band.GetNoDataValue()
if NDV is not None:
    # Convert No Data Points to nans
    Array[Array == NDV] = np.nan

# Now I'm ready to save the new file, in the meantime I have 
# closed the original, so I reopen it to get the projection
# information...
NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)

# Import new data
tex = gdal.Open('/Users/fasih/Dropbox/MobileOrg/textureShading-Asia/sichuan/east-0.01-land-only.png')
texArray = tex.ReadAsArray()
#NewArray = ((texArray.copy().astype(np.float) - 2**15*0) / 2**8).astype(np.uint8)
#DataType = GDT_Byte
NewArray = ((texArray.copy().astype(np.float) - 2**15*0) /
        2**0).astype(np.uint16)
DataType = GDT_UInt16

print NewArray.min(), NewArray.max()

# Set up the GTiff driver
driver = gdal.GetDriverByName('GTiff')

# Now turn the array into a GTiff.
NewFileName = CreateGeoTiff('newout', NewArray, driver, NDV, 
                            xsize, ysize, GeoT, Projection, DataType)
