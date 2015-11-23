import numpy as np
import textureShading as tex
filename = 'Data/east_0.1.tif'

data = tex.filenameToTexture(filename)

if True:
    np.save('texpy.npy', data)
    print "Saved raw"

data[tex.filenameToLandMask(filename)] = 0

tex.dataToPng(data, 'modularized.png', [0.5, 99.9], [.55, .4]) 
# also try [0.25, 0.2] for postPercentileScale. Or [0.55, 0.4]. Or [1, 1].
print "Saved PNG"
