import numpy
import scipy
from scipy import signal

def design(N = 32, passbandWidth = 0.03):
    if N % 2 != 0: raise ValueError('N must be even')
    if N < 2: raise ValueError('N must be > 1')
    if not (passbandWidth > 0 and passbandWidth < 0.5):
        raise ValueError('Need 0 < passbandWidth < 0.5')
    bands = numpy.array([0., .25 - passbandWidth, .25 + passbandWidth, .5])
    h = signal.remez(N + 1, bands, [1, 0], [1, 1])
    h[abs(h) <= 1e-4] = 0.0
    return h

def designAndSave(fname, *args, **kwargs):
    h = design(*args, **kwargs)
    h.tofile(fname, ",", "%.30g")
