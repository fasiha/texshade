import numpy as np
import gdal, gdalconst

# See http://www.gdal.org/classVRTRasterBand.html#a155a54960c422941d3d5c2853c6c7cef
def hist(fname, nbuckets=1000, percentiles=[1.0, 99.0]):
  """
  Given a filename, finds approximate percentile values and provides the
  gdal_translate invocation required to create an 8-bit PNG.

  Works by evaluating a histogram of the original raster with a large number of
  buckets between the raster minimum and maximum, then estimating the
  probability mass and distribution functions before reporting the percentiles
  requested.

  N.B. This technique is very approximate and hasn't been checked for asymptotic
  convergence. Heck, it uses GDAL's `GetHistogram` function in approximate mode,
  so you're getting approximate percentiles using an approximated histogram.

  Optional arguments:
  - `percentiles`: list of percentiles, between 0 and 100 (inclusive).
  - `nbuckets`: the more buckets, the better percentile approximations you get.
  """
  src = gdal.Open(fname)
  band = src.GetRasterBand(1)

  # Use GDAL to find the min and max
  (lo, hi, avg, std) = band.GetStatistics(True, True)

  # Use GDAL to calculate a big histogram
  rawhist = band.GetHistogram(min=lo, max=hi, buckets=nbuckets)
  binEdges = np.linspace(lo, hi, nbuckets+1)

  # Probability mass function. Triangle-intergration of this should yield 1.0.
  pmf = rawhist / (np.sum(rawhist) * np.diff(binEdges[:2]))
  # Cumulative probability distribution. Starts at 0, ends at 1.0.
  distribution = np.cumsum(pmf) * np.diff(binEdges[:2])

  # Which histogram buckets are close to the percentiles requested?
  idxs = [np.sum(distribution < p / 100.0) for p in percentiles]
  # These:
  vals = [binEdges[i] for i in idxs]

  # Append 0 and 100% percentiles (min & max)
  percentiles = [0] + percentiles + [100]
  vals = [lo] + vals + [hi]

  # Print the percentile table
  print "percentile (out of 100%),value at percentile"
  for (p, v) in zip(percentiles, vals):
    print "%f,%f" % (p, v)

  # Print out gdal_translate command (what we came here for anyway)
  print "gdal_translate -ot Byte -of PNG -a_nodata 0 -scale {min} {max} 0 255 {fname} clipped-{fname}.png".format(min=vals[1], max=vals[-2], fname=fname)


  return (vals, percentiles)


# Invoke as: `python hist.py my-raster.tif`.
if __name__ == '__main__':
  import sys

  if len(sys.argv) > 1:
    hist(sys.argv[1])
  else:
    print "python hist.py INPUT-RASTER"

