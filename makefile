Data/east_0.01.tif: Data/SRTM_SE_250m.tif Data/SRTM_NE_250m.tif
	gdalwarp -tr 0.01 0.01 -r bilinear Data/SRTM_SE_250m.tif Data/se_0.01.tif
	gdalwarp -tr 0.01 0.01 -r bilinear Data/SRTM_NE_250m.tif Data/ne_0.01.tif
	gdal_merge.py -o Data/east_0.01-raw.tif Data/se_0.01.tif Data/ne_0.01.tif
	gdal_translate -co COMPRESS=LZW Data/east_0.01-raw.tif east_0.01.tif

Data/east_0.1.tif: Data/east_0.01.tif
	gdalwarp -tr 0.1 0.1 -r bilinear Data/east_0.01.tif Data/east_0.1.tif

#Data/SRTM_SE_250m.tif: Data/SRTM_SE_250m_TIF.rar
#	unrar x -o- Data/SRTM_SE_250m_TIF.rar Data/
#
#Data/SRTM_NE_250m.tif: Data/SRTM_NE_250m_TIF.rar
#	unrar x -o- Data/SRTM_NE_250m_TIF.rar Data/
#
