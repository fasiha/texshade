from __future__ import division
import numpy as np
import sympy

fmax = 0.5
aS, vS, rS = sympy.symbols('a v r')
sol = sympy.integrate(rS ** aS * rS * sympy.besselj(0, 2 * sympy.pi * vS * rS),
                      (rS, 0, fmax)) * 2 * sympy.pi

alpha = 0.5

class EvalSol(object):
    def __init__(self, alpha):
        self.alpha = alpha
    def __call__(self, Ns):
        return np.array(map(lambda N: sol.evalf(subs={aS: self.alpha, vS: N}),
            Ns), dtype=np.float)
evaler = EvalSol(alpha)

radii = np.arange(0, 128*1024+1e-3, 1/32.0)
# radii = np.arange(0, 2*1024+1e-3, 1/32.0)
# radii = np.arange(0, 0.9, 1/32.0)

from multiprocessing import Pool
pool = Pool()
hyper = pool.map(evaler, np.array_split(radii, pool._processes))

hyper = np.hstack(hyper)

hyper.tofile('hyper.bin')
radii.tofile('r.bin')
