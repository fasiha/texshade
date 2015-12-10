import numpy as np
import gdal, gdalconst

import textureShading as tex
import big

origSize = np.array([28810, 100810])
filterSize = np.array([512, 512])

elevTif = '/Users/ahmed.fasih/Documents/Personal/textureShading-Asia/Data/SRTM_NE_250m.tif'
elevBin = '/Users/ahmed.fasih/Documents/Personal/textureShading-Asia/Data/ne-land-single.bin'
texBin = 'ne.bin'
texTif = 'ne.tif'

elevation = big.memmapInput(elevTif, elevBin, 'r')
elevationMask = elevation<-30000
del elevation


t = np.memmap(texBin, dtype=np.float32, mode='c', 
        shape=tuple(origSize + filterSize - 1))
tsub = t[filterSize[0]/2:filterSize[0]/2+origSize[0], 
        filterSize[1]/2:filterSize[1]/2+origSize[1] ]

percentiles = [0.1, 99.9]
postPercentileScale = [1, 1]

# Get histogram & calculate limits
skip = 100
nbins = 1000
hist, binEdges = np.histogram(tsub[np.logical_not(elevationMask)][::skip],
        bins=nbins, density=True)
distribution = (np.cumsum(hist) * np.diff(binEdges[:2]))
loIdx = np.sum(distribution < percentiles[0] / 100)
hiIdx = np.sum(distribution < percentiles[1] / 100)
if loIdx == 0 or hiIdx == nbins:
    print("Histogram with {} bins insufficient to estimate percentiles".format(nbins))
rawLimits = np.array([binEdges[loIdx], binEdges[hiIdx]])
limits = rawLimits * postPercentileScale

# Remap to output mask
tint = tex.touint(tsub, *tuple(limits), dtype=np.uint8)

# Reinstate mask
newNDV = 0
tint[elevationMask] = newNDV

# Save
driver = gdal.GetDriverByName('GTiff')
NDV, xsize, ysize, GeoT, Projection, DataType = tex.GetGeoInfo(elevTif)

tex.CreateGeoTiff(texTif, tint, driver, newNDV, xsize, ysize, GeoT, Projection,
        gdalconst.GDT_Byte)


"""
Then:

$ gdal_translate -of VRT -a_srs EPSG:4326 -gcp 0 0 -30.0104168  60.0104165 -gcp 100810 0 180.0104162 60.0104165 -gcp 100810 28810 180.0104162  -0.0104168 ne.tif ne.vrt
$ gdalwarp -of VRT -t_srs EPSG:4326 ne.vrt ne2.vrt
$ gdal2tiles.py -p mercator -a 0,0,0 ne2.vrt

But the basemap doesn't seem to agree. And the ne2.vrt doesn't load in QGIS.
"""

