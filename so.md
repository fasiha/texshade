# Background
In image processing applications like [elevation texture shading](http://textureshading.com/Home.html), one applies a sharpening filter that is defined in the frequency domain as:
$$
J = \mathscr{F}^{-1} \left[ \left| f \right|^\alpha \cdot \mathscr{F}\left[I\right] \right],
$$
where $0 \leq \alpha \leq 2$. That is, 

1. take the 2D FFT $\mathscr{F}$ of an image $I$, 
2. weight each Fourier bin by the fractional root of the radial spatial frequency $\left| f \right|^\alpha = \left( f_x^2 + f_y^2 \right)^{\frac{\alpha}{2}}$, and
3. take the inverse 2D FFT $\mathscr{F}^{-1}$.

In Matlab/Octave,

    alpha = 0.5;
    % load I
    [Ny, Nx] = size(I); % assume Nx and Ny are even
    fx = [-Nx / 2 : Nx / 2 - 1] / Nx;
    fy = [-Ny / 2 : Ny / 2 - 1]' / Ny;
    J = ifft2(ifftshift(bsxfun(@plus, fy.^2, fx.^2).^(alpha / 2) .* fftshift(fft2(I))))

(For the Python/R users: `bsxfun` performs array broadcasting, in this case adding each element of a row vector `fx.^2` to each element of a column vector `fy.^2` to make a new $N_y \times N_x$ array. Equivalently, I could have something like `[fxMat, fyMat] = meshgrid(fx, fy);` and replaced the call to `bsxfun` with something like `(fxMat.^2 + fyMat.^2)`.)

# Problem
I want to run this algorithm on extremely large images `I`, larger than computer memory, e.g., with memory-mapped files. This precludes using a 2D FFT, because standard FFT libraries would need enough memory to store their outputs.

# Approach
We know that a filter with a certain frequency response has a spatial-domain impulse response, and convolving this impulse response with the image in the spatial domain can be equivalent to frequency-domain filtering.

To work around the memory constraints, I'd like to approximate the filtering operation by truncating the impulse response: instead of a 100k by 100k impulse response for an image of the same size, truncate to maybe a 1000 by 1000 filter, and apply overlap-add filtering on small chunks of the input image—and generating the output in small chunks as well.

The simplest way to find the impulse response, 2D IFFT of $|f|^\alpha$, i.e., `ifft2(ifftshift(bsxfun(@plus, fy.^2, fx.^2).^(alpha / 2)))`, is obviously not going to work for large filters for the same reasons as above. So I'm looking for a way to find the spatial-domain impulse response of this fractional-root frequency-response filter.

# Analysis

We know that a circularly symmetric function, i.e., one that is in polar coordinates depends only on the radius: $g(x,y) = g\left(\sqrt{x^2 + y^2}\right)$, also has a circularly symmetric discrete Fourier transform (DFT). Furthermore, this circularly symmetric DFT is a [Hankel transform of the zeroth order](https://en.wikipedia.org/wiki/Hankel_transform#Relation_to_the_Fourier_transform_.28circularly_symmetric_case.29).

Unlike Fourier transforms, I'm at all familiar with Hankel transforms, beyond what I've studied in the last couple of days. But assuming a 1D Hankel transform is computable with FFT-like methods, I could obtain the impulse response of the fractional-root filter $|f|^\alpha$ in one dimension, along spatial radius, without onerous memory requirements.

# What I've tried

I downloaded Adam Wyatt's [Hankel Transform](http://www.mathworks.com/matlabcentral/fileexchange/15623-hankel-transform) Matlab package. It implements this [paper by Guizar-Sicairos & Gutierrez-Vega](http://www.acms.arizona.edu/FemtoTheory/MK_personal/opti547/literature/QD-Hankel.pdf). This implementation of the discrete Hankel transform (DHT) won't meet my needs because it evaluates an N-point Hankel transform by first generating an N by N matrix, but I'd like to verify that the principle works.

Download the above code and **put the code in a directory called “+Hankel” in your path** (my code below uses namespaces, and in Matlab, functions inside directories prepended with a plus sign are namespaces).

I can reproduce Figure 1(a) and 1(c) from Guizar-Sicairos & Gutierrez-Vega's paper above using Wyatt's Hankel transform code. This figure shows the fourth order, 256-point discrete Hankel transform of a sinc function defined on $[0, 3]$:

    pOrder = 4; % change to 1 to get Figure 1(a)
    R = 3;
    Nr = 256;
    h = Hankel.hankel_matrix(pOrder, R, Nr);

    % parameter for the test function, a sinc
    gamma = 5;
    % Our sinc function includes pi: sinc(x) = x == 0 ? 1 : sin(pi * x) / (pi * x)
    f1fun = @(r) sinc(2 * gamma * r);

    % Samples of f1(r) are on the zeros of Bessel function
    r = h.J_roots / (2 * pi * h.vmax);

    % Intermediate forms
    F1 = f1fun(r) ./ h.J * h.vmax;
    % Applying the linear transform operator
    F2 = h.T * F1;

    % Transformed vector
    f2 = F2 .* h.J / h.vmax;
    % Transformed vector's sample points
    v = h.J_roots / (2 * pi * h.rmax);

    figure; plot(v, f2)

([Gist](https://gist.github.com/fasiha/16f796be6a55cc7dae98) if you prefer cloning, but that doesn't use the `Hankel` namespace, since gists don't support directories.)

[![Figure 1(c)][1]][1]

And in fact, if you change `pOrder = 0`, we get the circularly symmetric DFT along radius: a sinc at frequency of `gamma = 5` Hz in the spatial domain corresponds to a rectangular window of width 5 in the frequency domain (with the characteristic spike at radial frequency of 5 caused by fininte truncation of the sinc):

[![Using 0th order Hankel transform, i.e., the circularly symmetric DFT][2]][2]

From these examples involving sincs, I can tell the DHT is working. But when I try changing the example above to use `f1fun = @(r) r.^0.5;`, i.e., evaluate $r^\alpha = r^{0.5}$, I get nonsensical results: change the beginning of the code snippet above to:

    pOrder = 0;
    R = 1;
    Nr = 256;
    h = Hankel.hankel_matrix(pOrder, R, Nr);
    alpha = 0.5;
    f1fun = @(r) r.^alpha;
    % ... rest as above ...

and reusing the rest of the code, produces this Hankel transform:

***Mostly solved*** **SEE sympystuff.py** for a working model. However, it needs to enforce the full |f| less than R_radius, otherwise it's not happy. That means the image will be somewhat LPF'd. See **tex.m** for visual impact of this radial LPF. It's not too bad actually. But you can clearly see the ringing introduced by it. The interpolation filters that Leland Brown suggests (approximating the cubic interpolator in the frequency domain) do help reduce the ringing, without any painful blurring (like Hamming windowing the frequency data did). So it should be possible to use the hypergeometric expression found in `sympystuff.py` as well as a cubic interpolating spatial-domain filter, combine the two into a single FIR filter (truncated appropriately). Oh, note that the cubic spline interpolator is separable.

Leland Brown also suggests downsampling the big image, texture-shading that, and somehow using that to guide stitching together of smaller texture-shaded images, but this simply papers-over the problems at the edges, and it needs you to filter the big image to downsample it.

  [1]: http://i.stack.imgur.com/kL9lr.png
  [2]: http://i.stack.imgur.com/QtRxr.png
