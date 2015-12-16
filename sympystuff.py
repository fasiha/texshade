from __future__ import division
import numpy as np
import sympy
import scipy.interpolate as interpolate

import pylab
pylab.ion()


def mkim(x):
    pylab.figure()
    pylab.imshow(x)
    pylab.colorbar()
    pylab.grid()


def vec(x):
    return x.reshape([x.size, -1])


fmax = 0.5
aS, vS, rS = sympy.symbols('a v r')
sol = sympy.integrate(rS ** aS * rS * sympy.besselj(0, 2 * sympy.pi * vS * rS),
                      (rS, 0, fmax)) * 2 * sympy.pi

alpha = 0.5
N = 20

if True:
    radii = np.arange(0, 32*1024+1e-3, 1/32.0)
    from multiprocessing import Pool
    pool = Pool(processes=7)
    hyper = pool.map(lambda x: sol.evalf(subs={aS: alpha, vS: x}), radii)

evalSol = np.vectorize(lambda x: sol.evalf(subs={aS: alpha, vS: x}),
        otypes=[float])
denseVec = np.linspace(0, N / 2, 256)
solDense = evalSol(denseVec)

import numpy.fft as fft
mkf = lambda N: np.arange(-N / 2.0, N / 2.0) / N
mkt = lambda N: np.arange(-N / 2.0, N / 2.0)
x = vec(mkt(N)).T
y = vec(mkt(N))
r = np.sqrt(x ** 2 + y ** 2)
spatialFilter = evalSol(r)

freqResp = np.real(fft.fftshift(fft.fft2(fft.ifftshift(spatialFilter))))
freqRespRow = fft.fftshift(fft.ifftshift(freqResp)[0])

fx = vec(mkf(N)).T
fy = vec(mkf(N))
fr = np.sqrt(fx ** 2 + fy ** 2)

pylab.figure()
pylab.plot(fx.ravel(), freqRespRow)
pylab.title('frequency domain: obviously wrong, not |f|^alpha')

freqRespCalc = (fr ** alpha) * (fr <= fmax)
spatialFilterCalc = np.real(fft.fftshift(fft.fft2(
    fft.ifftshift(freqRespCalc)))) / freqRespCalc.size
spatialCalcRow = fft.fftshift(fft.ifftshift(spatialFilterCalc)[0])

remmax = lambda x: x / np.max(x)
pylab.figure()
pylab.plot(x.ravel(), spatialCalcRow, 'o-', denseVec, solDense, x.ravel(),
           evalSol(abs(x.ravel())), 'og')
pylab.legend('calculated,theoretical,discrete theoretical'.split(','))
pylab.grid()
pylab.title('spatial domain')
