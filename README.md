## Setup
Download a bunch of ZIP files.

Uncompress them all.

Either:
- Make a copy of each one, burning NDV, lakes, and oceans to 0.
- Make one big ENVI.
- Pros: parallelization. Cons: 400 GB more needed.
OR
- Make one big VRT,
- and make one big float32 ENVI burning NDV, lakes, and oceans to 0.

## Texture generation
Run through texture generator & save the texture as a new binary file of float32.

## Post-processing
Burn NDVs in the texture bin for anywhere the the original input was NDV or lake or ocean. (Two passes: find NDVs in the original elevation, and then lakes and coastlines.)

Find 1th and 99th percentiles (or something similar) for the non-NDV texture values.

Map the float32 texture bin to byte, interpolated between the lower and upper thresholds. Let NDV be NDV.

Make tiles using the byte texture.

