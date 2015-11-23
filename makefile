Data/east_0.01.tif: Data/SRTM_SE_250m.tif Data/SRTM_NE_250m.tif Data/ne_10m_land/ne_10m_land.shp
	gdalwarp -tr 0.01 0.01 -r average Data/SRTM_SE_250m.tif Data/se_0.01.tif
	gdalwarp -tr 0.01 0.01 -r average Data/SRTM_NE_250m.tif Data/ne_0.01.tif
	gdal_merge.py -n -32768 -a_nodata -32768 -o Data/east_0.01.tif Data/se_0.01.tif Data/ne_0.01.tif
	# Crop to land
	gdalwarp -cutline Data/ne_10m_land/ne_10m_land.shp -dstnodata -32768 Data/east_0.01.tif Data/land.tif
	# Compress
	gdal_translate -a_nodata -32768 -co COMPRESS=LZW Data/land.tif Data/east_0.01.tif

Data/east_0.1.tif: Data/east_0.01.tif Data/ne_10m_land/ne_10m_land.shp
	# Downsample the downsampled image
	gdalwarp -tr 0.1 0.1 -r average Data/east_0.01.tif Data/east_0.1.tif
	# Crop to land
	gdalwarp -cutline Data/ne_10m_land/ne_10m_land.shp -dstnodata -32768 Data/east_0.1.tif Data/land.tif
	# Rename
	mv Data/land.tif Data/east_0.1.tif

Data/ne_10m_land/ne_10m_land.shp:
	curl -O http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_land.zip
	mkdir -p Data/ne_10m_land
	mv ne_10m_land.zip Data/ne_10m_land
	cd Data/ne_10m_land
	unzip ne_10m_land.zip
	cd ../..

#Data/SRTM_SE_250m.tif: Data/SRTM_SE_250m_TIF.rar
#	unrar x -o- Data/SRTM_SE_250m_TIF.rar Data/
#
#Data/SRTM_NE_250m.tif: Data/SRTM_NE_250m_TIF.rar
#	unrar x -o- Data/SRTM_NE_250m_TIF.rar Data/
#
