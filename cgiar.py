import numpy as np
import textureShading as tex
filename = 'Data/east_0.01.tif'

generateData = True # if false, tries to load raw from disk
saveRaw = True
writeTif = True
writePng = True

if generateData:
    data = tex.filenameToTexture(filename)
else:
    data = np.load('texpy.npy')
    print "Loaded data"

if saveRaw and not generateData: # if data wasn't generated, it was loaded
    np.save('texpy.npy', data)
    print "Saved raw"

# Re-apply mask
data[tex.filenameToLandMask(filename)] = 0

if writeTif:
    tex.dataToGTiff16(data, filename, 'modularized.tif', [0.5, 99.9], [0.55, 0.4])
    print "Saved TIF"

if writePng:
    tex.dataToPng16(data, 'modularized.png', [0.5, 99.9], [.55, .4]) 
    # also try [0.25, 0.2] for postPercentileScale. Or [0.55, 0.4]. Or [1, 1].
    print "Saved PNG"

