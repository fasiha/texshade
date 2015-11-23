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



# filename = '/Users/fasih/Dropbox/Data/cgiar/SRTM_NE_250m.tif'
filename = '/Users/fasih/Dropbox/Data/cgiar/0.01_small_East.tif'
# filename = '/Users/fasih/Dropbox/Data/chgis/chgis_dem.tif'
fileHandle = gdal.Open(filename, gdalconst.GA_ReadOnly)
band = fileHandle.GetRasterBand(1)
data = band.ReadAsArray().astype(np.float32)
if not False:
    data = np.maximum(0, data)
print "Read data"

data = sciMoreReal(data, 0.5)
print "Computed texture"

np.save('texpy.npy', data)
print "Saved and exiting"
