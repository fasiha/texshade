from multiprocessing import cpu_count
NCPUS = cpu_count()

import numpy as np
from numpy.fft import fft2, ifft2, rfft2, irfft2
import pyfftw
pyfftw.interfaces.cache.enable()


def overlapadd2(A, Hmat, L=None, Nfft=None, y=None, verbose=False, Na=None):
    """
    Fast two-dimensional linear convolution via the overlap-add method.

    The overlap-add method is well-suited to convolving a very large array,
    `A`, with a much smaller filter array, `Hmat` by breaking the large
    convolution into many smaller `L`-sized sub-convolutions, and evaluating
    these using the FFT. The computational savings over the straightforward
    two-dimensional convolution via, say, scipy.signal.convolve2d, can be
    substantial for large A and/or Hmat.

    Parameters
    ----------
    A, Hmat : array_like
        Two-dimensional input arrays to be convolved. For computational
        purposes, A should be larger.
    L : sequence of two ints, optional
        Length of sub-convolution to use. This should ideally be as large as
        possible to obtain maximum speedup: that is, if you have the memory to
        compute the linear convolution using FFT2, i.e., replacing spatial
        convolution with frequency-domain multiplication, then let `L =
        np.array(A.shape) + np.array(Hmat.shape) - 1`. Usually, though, you
        are considering overlap-add because you can't afford a batch transform
        on your arrays, so make `L` as big as you can.
    Nfft : sequence of two ints, optional
        Size of the two-dimensional FFT to use for each `L`-sized
        sub-convolution. If omitted, defaults to the next power-of-two, that
        is, the next power-of-two on or after `L + np.array(Hmat.shape) - 1`.
        If you choose this default, try to avoid `L` such that this minimum, `L
        + np.array(Hmat.shape) - 1`, is just a bit greater than a power-of-two,
        since you'll be paying for a 4x larger FFT than you need.
    y : array_like, optional
        Storage for the output. Useful if using a memory-mapped file, e.g.
    verbose : boolean, optional
        If True, prints a message for each `L`-sized subconvolution.

    Returns
    -------
    y : same as passed in, or ndarray if no `y` passed in
        The `np.array(A.shape) + np.array(Hmat.shape) - 1`-sized
        two-dimensional array containing the linear convolution. Should be
        within machine precision of, e.g., `scipy.signal.convolve2d(A,
        Hmat, 'full')`.

    Raises
    ------
    ValueError if `L` and `Nfft` aren't two-element, and too small: both
    elements of `L` must be greater than zero, and `Nfft`'s must be greater
    than `L + np.array(Hmat.shape) - 1`. Also if `A` or `Hmat` aren't
    two-dimensional arrays, or if `y` doesn't have the correct size to store
    the output of the linear convolution.

    References
    ----------
    Wikipedia is only semi-unhelpful on this topic: see "Overlap-add method".
    """
    M = np.array(Hmat.shape)
    if not hasattr(A, '__call__'): # not a function!
        Afun = lambda a, b, c, d: A[a:b, c:d]
    else:
        Afun = A

    if Na is None:
        Na = np.array(A.shape)

    yIsFunction = False
    if y is None:
        y = np.zeros(M + Na - 1, dtype=A.dtype)
    elif hasattr(y, '__call__'):
        yIsFunction = True
    elif y.shape != tuple(M + Na - 1):
        raise ValueError('y given has incorrect dimensions', M + Na - 1)

    if L is None:
        L = M * 100
    else:
        L = np.array(L)

    if Nfft is None:
        Nfft = 2 ** np.ceil(np.log2(L + M - 1)).astype(int)
    else:
        Nfft = np.array(Nfft, dtype=int)

    if not (np.all(L > 0) and L.size == 2):
        raise ValueError('L must have two positive elements')
    if not (np.all(Nfft >= L + M - 1) and Nfft.size == 2):
        raise ValueError('Nfft must have two elements >= L + M - 1 where '
                         'M = Hmat.shape')

    realOnly = np.isrealobj(Afun(0, 1, 0, 1)) and np.isrealobj(Hmat)
    if realOnly:
        xAligned = pyfftw.zeros_aligned(tuple(Nfft), dtype='float32')
        (myifft2, myfft2) = (pyfftw.interfaces.numpy_fft.irfft2, 
                pyfftw.interfaces.numpy_fft.rfft2)
    else:
        assert(False)

    HmatAligned = pyfftw.empty_aligned(tuple(Nfft), dtype='float32')
    HmatAligned[:Hmat.shape[0], :Hmat.shape[1]] = Hmat
    Hf = pyfftw.interfaces.numpy_fft.rfft2(Hmat, Nfft)
    del HmatAligned

    (XDIM, YDIM) = (1, 0)
    start = [0, 0]
    endd = [0, 0]
    while start[XDIM] <= Na[XDIM]:
        endd[XDIM] = min(start[XDIM] + L[XDIM], Na[XDIM])
        start[YDIM] = 0
        while start[YDIM] <= Na[YDIM]:
            endd[YDIM] = min(start[YDIM] + L[YDIM], Na[YDIM])
            if verbose:
                print("Input start", start, "Input end", endd, "FFT size", Nfft)
            
            diffs = [endd[YDIM] - start[YDIM], endd[XDIM] - start[XDIM]]
            xAligned[:] = 0
            xAligned[:diffs[0], :diffs[1]] = Afun(start[YDIM], endd[YDIM], start[XDIM], endd[XDIM])
            yt = myifft2(Hf * myfft2(xAligned, threads=NCPUS/2), threads=NCPUS/2)

            thisend = np.minimum(Na + M - 1, start + Nfft)
            if yIsFunction:
                y(start[YDIM], thisend[YDIM], start[XDIM], thisend[XDIM],
                        (yt[:(thisend[YDIM] - start[YDIM]), :(thisend[XDIM] - start[XDIM])]))
            else:
                y[start[YDIM] : thisend[YDIM], start[XDIM] : thisend[XDIM]] += (
                    yt[:(thisend[YDIM] - start[YDIM]),
                        :(thisend[XDIM] - start[XDIM])])

            start[YDIM] += L[YDIM]

            del yt  # helps!
        start[XDIM] += L[XDIM]
    return y


def test():
    from scipy.signal import convolve2d
    A = np.random.randn(33, 55)
    H = np.random.randn(4, 5)
    gold = convolve2d(A, H)
    assert(np.allclose(gold, overlapadd2(A, H, L=[12, 12])))
    assert(np.allclose(gold, overlapadd2(A, H, L=[12, 120])))
    assert(np.allclose(gold, overlapadd2(A, H, L=[90, 120])))

    assert(np.allclose(gold, overlapadd2(H, A, L=[190, 220])))
    assert(np.allclose(gold, overlapadd2(H, A, L=[1, 1])))

    assert(np.allclose(gold, overlapadd2(H, A)))
    assert(np.allclose(gold, overlapadd2(A, H)))

    assert(np.allclose(convolve2d(A.T, H.T),
                       overlapadd2(A.T, H.T, L=[190, 220])))

    A = np.random.randn(33, 55) + 1j * np.random.randn(33, 55)
    H = np.random.randn(4, 5) + 1j * np.random.randn(4, 5)
    gold = convolve2d(A, H)
    assert(np.allclose(gold, overlapadd2(H, A)))
    assert(np.allclose(gold, overlapadd2(A, H)))

    print('Test passed!')


if __name__ == '__main__':
    test()


