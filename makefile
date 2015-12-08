Data/ne-land-single.bin: Data/SRTM_NE_250m.tif Data/ne_10m_land/ne_10m_land.shp
	gdalwarp -cutline Data/ne_10m_land/ne_10m_land.shp -dstnodata -32768 Data/SRTM_NE_250m.tif Data/ne-land.tif
	gdal_translate -of ENVI -ot Float32 Data/ne-land.tif Data/ne-land-single.bin
	rm Data/ne-land.tif

Data/koeppen-0.1.tif Data/koeppen-0.01.tif: Data/Koeppen-Geiger-GIS/koeppen-geiger.shp
	gdal_rasterize -co COMPRESS=LZW -a "GRIDCODE"  -a_nodata 0 -tr 0.0100000000000000 0.0100000000000000 -ot Byte -te -30.0104168 -60.0095835 180.0095832 60.0104165 Data/Koeppen-Geiger-GIS/koeppen-geiger.shp Data/koeppen-0.01.tif
	gdal_rasterize -co COMPRESS=LZW -a "GRIDCODE"  -a_nodata 0 -tr 0.100000000000000 0.100000000000000 -ot Byte -te -30.0104168 -59.9895835 179.9895832 60.0104165  Data/Koeppen-Geiger-GIS/koeppen-geiger.shp Data/koeppen-0.1.tif

Data/east-land.tif Data/west-land.tif: Data/SRTM_SE_250m.tif Data/SRTM_NE_250m.tif Data/ne_10m_land/ne_10m_land.shp Data/SRTM_W_250m.tif
	# not working, make dies at 70%: # gdal_merge.py -n -32768 -a_nodata -32768 -co COMPRESS=LZW -o Data/all.tif Data/SRTM_SE_250m.tif Data/SRTM_NE_250m.tif	Data/SRTM_W_250m.tif
	# Crop to land: don't bother with compression
	gdalwarp -cutline Data/ne_10m_land/ne_10m_land.shp -dstnodata -32768 Data/SRTM_SE_250m.tif Data/se-land.tif
	gdalwarp -cutline Data/ne_10m_land/ne_10m_land.shp -dstnodata -32768 Data/SRTM_NE_250m.tif Data/ne-land.tif
	gdalwarp -cutline Data/ne_10m_land/ne_10m_land.shp -dstnodata -32768 Data/SRTM_W_250m.tif Data/w-land.tif
	# Merge NE and SE into one so they can be processed together.
	gdal_merge.py -n -32768 -a_nodata -32768 -o Data/east-land.tif Data/ne-land.tif Data/se-land.tif
	# delete NE and SE lands
	rm -f Data/ne-land.tif
	rm -f Data/se-land.tif


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

Data/ne2-0.1.tif: Data/NE2_HR_LC/NE2_HR_LC.tif
	gdalwarp -tr 0.1 0.1 -r average Data/NE2_HR_LC/NE2_HR_LC.tif Data/ne2-0.1.tif
	gdaltindex shp-0.1.shp Data/east_0.1.tif
	gdalwarp -dstnodata 0 -q -cutline shp-0.1.shp -crop_to_cutline -of GTiff Data/ne2-0.1.tif ne2-0.1-east.tif

Data/ne2-0.01.tif: Data/NE2_HR_LC/NE2_HR_LC.tif
	gdalwarp -tr 0.01 0.01 -r average Data/NE2_HR_LC/NE2_HR_LC.tif Data/ne2-0.01.tif
	gdaltindex shp-0.01.shp Data/east_0.01.tif
	gdalwarp -q -cutline shp-0.01.shp -crop_to_cutline -of GTiff Data/ne2-0.01.tif ne2-0.01-east.tif
