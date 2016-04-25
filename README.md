## Setup
Download a bunch of ZIP files.

Uncompress them all.

Make one big VRT.

## Texture generation
(Make the Hankel look-up table. Took 44 minutes for 0:1/32.0:128K and 32 cores.)

Run through texture generator & save the texture as a new binary file of float32.
- real    1843m50.867s, user    416m54.831s, sys     267m26.899s
- real    1699m42.307s, user    519m30.430s, sys     276m52.106s
- real    1934m58.130s, user    547m27.940s, sys     342m29.861s
- not sure which is which, but the tex5-clean-pyfftw.bin finished first, despite being started much later. But all three finished within 10 minutes of each other.

## Post-processing
Burn NDVs in the texture bin for anywhere the the original input was NDV or lake or ocean. (Two passes: find NDVs in the original elevation, and then lakes and coastlines.)

Find 1th and 99th percentiles (or something similar) for the non-NDV texture values.

Map the float32 texture bin to byte, interpolated between the lower and upper thresholds. Let NDV be NDV.

Make tiles using the byte texture.

