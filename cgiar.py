import numpy as np
import textureShading as tex
filename = 'Data/east_0.01.tif'
filename = 'Data/east_0.1.tif'

generateData = True # if false, tries to load raw from disk
saveRaw = True
writeTif = True
writePng = not True

if generateData:
    data = tex.filenameToTexture(filename)
else:
    data = np.load('texpy.npy')
    print "Loaded data"

if saveRaw and not generateData: # if data wasn't generated, it was loaded
    np.save('texpy.npy', data)
    print "Saved raw"

# Re-apply mask
data[tex.filenameToNoDataMask(filename)] = 0

percentiles = [0.5, 99.9]
postPercentileScale = [1, 1.0]
postPercentileScale = [0.25, 0.2]
postPercentileScale = [0.55, 0.4]
if writeTif:
    tex.dataToGTiff16(data, filename, 'modularized.tif', percentiles, postPercentileScale)
    print "Saved TIF"

if writePng:
    tex.dataToPng16(data, 'modularized.png', percentiles, postPercentileScale) 
    # also try [0.25, 0.2] for postPercentileScale. Or [0.55, 0.4]. Or [1, 1].
    print "Saved PNG"

