#!/N/u/macblab/Carbonate/.conda/envs/py37_env/bin/python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import glob, math
from osgeo import gdal
import time
start = time.time()

#Getting list of reflectance files 
files_to_mosaic = glob.glob('/N/project/SRER/Unmixing/BRDFHeight/10TSig/mask/FC/*EM2.tif')
#Defining a function to get extent of each tiff file
def get_extent(fn):
    "'returns min_x, max_y, max_x, min_y'"
    ds = gdal.Open(fn)
    gt = ds.GetGeoTransform()
    return (gt[0], gt[3], gt[0]+gt[1]*ds.RasterXSize, gt[3]+gt[5]*ds.RasterYSize)

#getting the extent of 1st file in the classified file list 
min_x, max_y, max_x, min_y = get_extent(files_to_mosaic[0])

#Getting the extent of mosaic file
for fn in files_to_mosaic[1:]:
    minx, maxy, maxx, miny = get_extent(fn)
    min_x = min(min_x, minx)
    max_y = max(max_y, maxy)
    max_x = max(max_x, maxx)
    min_y = min(min_y, miny)

#Getting the number of rows and columns needed to creat mosaic file
in_ds = gdal.Open(files_to_mosaic[0])
gt = in_ds.GetGeoTransform()
rows = math.ceil((max_y-min_y)/-gt[5])
columns = math.ceil((max_x-min_x)/gt[1])
NP2GDAL_CONVERSION = {
        "uint8": 1,
        "int8": 1,
        "uint16": 2,
        "int16": 3,
        "uint32": 4,
        "int32": 5,
        "float32": 6,
        "float64": 7,
        "complex64": 10,
        "complex128": 11,}
gdaltype = NP2GDAL_CONVERSION[in_ds.GetRasterBand(1).ReadAsArray().dtype.name]
#Creating a blank mosaic file
path= '/N/project/SRER/Merged/FC/'
file_name = path + 'Supervised_EM2_Grass.tif'
driver = gdal.GetDriverByName("GTiff")
out_ds = driver.Create(file_name, columns, rows, 1, gdaltype)
out_ds.SetProjection(in_ds.GetProjection())
#out_band = out_ds.GetRasterBand(1)

#georeferencing the mosaic file
gt = list (in_ds.GetGeoTransform())
gt[0], gt[3]=min_x, max_y
out_ds.SetGeoTransform(gt)
"""reading each reflectance file as numpy nd array. On the same location of that file inside the mosaic file there arrays are written according to their bands"""
for fn in files_to_mosaic:
    #getting the location where to write the data
    in_ds = gdal.Open(fn)
    gt = in_ds.GetGeoTransform()
    minx=gt[0]
    maxy=gt[3]
    x=math.ceil(minx-min_x)
    y=math.ceil(-(maxy-max_y))
    #writing each array in each band
    out_ds.GetRasterBand(1).WriteArray(in_ds.GetRasterBand(1).ReadAsArray(),x,y)
    print(fn)
del in_ds, out_ds

end = time.time()
print("tot time: ", (end-start))
