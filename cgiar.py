import numpy as np
import textureShading as tex
filename = 'Data/east_0.1.tif'
filename = 'Data/east_0.01.tif'

generateData = True # if false, tries to load raw from disk
saveRaw = True
writeTif = True
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

# NE2 blending
ne2file = 'ne2-0.1-east.tif'
# ne = tex.filenameToData(ne2file, None)
