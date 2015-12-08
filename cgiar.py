import numpy as np
import textureShading as tex

mksmall = lambda : ('Data/east_0.1.tif', 'Data/koeppen-0.1.tif')
mkmed = lambda : ('Data/east_0.01.tif', 'Data/koeppen-0.01.tif')

filename, koeppenRasterFile = mksmall()

generateData = not True # if false, tries to load raw from disk
saveRaw = True
writeTif = not True
writePng = not True

if generateData:
    data = tex.filenameToTexture(filename)
else:
    data = np.load('texpy.npy')
    print "Loaded saved data"

if saveRaw and generateData: # if data wasn't generated, it was loaded
    np.save('texpy.npy', data)
    print "Saved raw"

# Re-apply mask
data[tex.filenameToNoDataMask(filename)] = 0

# Some presets for making pretty pictures right off the bat
percentiles = [0.5, 99.9]
postPercentileScale = [1, 1.0]
postPercentileScale = [0.25, 0.2]
postPercentileScale = [0.55, 0.4]

# A preset for images whose curves we'll later manipulate: so save a lot of data
percentiles = [0.1, 99.9]
postPercentileScale = [1, 1]

if writeTif:
    tex.dataToGTiff(data, filename, 'modularized.tif', 16/2, percentiles,
            postPercentileScale)
    print "Saved TIF"

if writePng:
    tex.dataToPng16(data, 'modularized.png', percentiles, postPercentileScale) 
    print "Saved PNG"

# Koeppen-Geiger climate coloring
if 1:
    # This is the min/max of the texture we'll be using
    dataMask = np.logical_not(tex.filenameToNoDataMask(filename))
    limits = tex.dataToPercentileLimits(data[dataMask], 
                                        percentiles, postPercentileScale)
    # Here's the 8-bit texture: 0 to 255. These will be remapped.
    tex8bit = tex.touint(data, *limits, dtype=np.uint8)
    tex8bit[np.logical_not(dataMask)] = 0

    # Here are the perceptually-balanced colors, five channels for five climate
    # zones
    import csv
    import matplotlib.colors as mplcolors
    koeppenCsvFile = 'Data/koeppen-colors-low-contrast.csv'
    with open(koeppenCsvFile, 'r') as csvfile:
        reader = csv.reader(csvfile)
        colors = [np.array(map(mplcolors.hex2color, r)) for r in reader]
    legend = "A,D,E,C,B".split(',')
    legendDict = dict(zip(legend, range(len(legend))))
    numColors = len(colors[0])

    def intZoneToColor(values, zone):
        """
        values: uint8 array_like, zone: 'A' through 'E'
        """
        return colors[legendDict[zone]][((np.array(values) * (numColors - 1e-6))
                                          / 255.0).astype(int)]
    
    # Here's how to handle the Koeppen-Geiger raster we've created (from the
    # original SHP file): convert its values to indexes into `colors`
    koeppenShpNumbers = map(int,
            '11,12,13,14,21,22,26,27,31,32,33,34,35,36,37,38,39,41,42,43,44,45,46,47,48,49,50,51,52,61,62'.split(','))
    koeppenShpCodes = 'Af,Am,As,Aw,BWk,BWh,BSk,BSh,Cfa,Cfb,Cfc,Csa,Csb,Csc,Cwa,Cwb,Cwc,Dfa,Dfb,Dfc,Dfd,Dsa,Dsb,Dsc,Dsd,Dwa,Dwb,Dwc,Dwd,EF,ET'.split(',')
    koeppenShpToLegend = dict(zip(koeppenShpNumbers, 
                                  [legendDict[c[0]] for c in koeppenShpCodes]))

    # Load the Koeppen-Geiger raster
    koeppenRaster = tex.filenameToData(koeppenRasterFile, np.int8)

    # Make the masks for each of the five climates
    intervalInclusive = lambda x, lo, hi: np.logical_and(x >= lo, x <= hi)
    koeppenA = intervalInclusive(koeppenRaster, 11, 14)
    koeppenB = intervalInclusive(koeppenRaster, 21, 27)
    koeppenC = intervalInclusive(koeppenRaster, 31, 39)
    koeppenD = intervalInclusive(koeppenRaster, 41, 52)
    koeppenE = intervalInclusive(koeppenRaster, 61, 62)
    
    dataClim = np.dstack([tex8bit, tex8bit, tex8bit])
    # output of intZoneToColor is a float between [0, 1]. To output this as a
    # byte, scale it by this (256 - epsilon), so [0, 1] -> [0, 255] inside the
    # byte-array. We do 256 - epsilon so that 1 gets scaled to nearly 256, which
    # becomes 255 when floored by the float-to-int. This way, all 256 levels of
    # dynamic range are used. We don't want to scale this by 255 because only
    # pixel values of 1.0 will be mapped to 255, instead of pixels with values
    # very close to 1.0.
    oneTo255 = (256 - 1e-6)
    dataClim[koeppenA, :] = intZoneToColor(tex8bit[koeppenA], 'A') * oneTo255
    dataClim[koeppenB, :] = intZoneToColor(tex8bit[koeppenB], 'B') * oneTo255
    dataClim[koeppenC, :] = intZoneToColor(tex8bit[koeppenC], 'C') * oneTo255
    dataClim[koeppenD, :] = intZoneToColor(tex8bit[koeppenD], 'D') * oneTo255
    dataClim[koeppenE, :] = intZoneToColor(tex8bit[koeppenE], 'E') * oneTo255
    
    dataClim[np.logical_not(dataMask), :] = 0

    # Write GeoTIFF
    import gdal, gdalconst
    driver = gdal.GetDriverByName('GTiff')
    NDV, xsize, ysize, GeoT, Projection, DataType = tex.GetGeoInfo(filename)
    newNDV = 0
    tex.CreateGeoTiff('col.tif', dataClim, driver, newNDV, xsize, ysize, 
            GeoT, Projection, gdalconst.GDT_Byte, numBands=3)

