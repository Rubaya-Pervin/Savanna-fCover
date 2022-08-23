#!/N/u/macblab/Carbonate/.conda/envs/py37_env/bin/python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import glob, math, subprocess, h5py
from osgeo import gdal

#Getting list of CHM and reflectance files 
files_to_mosaic1 = glob.glob('N*CHM.tif')
files_to_mosaic2 = glob.glob('*.h5')

"""the number of CHM files are larger than the number of reflectance files covering larger area. For our purporse we need same extent of CHM and reflectance data. As the file names contain the coordinates of the files we are only taking the files that have same coordinates from all file list of CHM and reflectance""" 
#Creating a list of common CHM and Reflectance file
files_to_masaic_height = []
files_to_masaic_reflct = []
for height in files_to_mosaic1:
    sep_h =height.split('_')
    for reflct in files_to_mosaic2:
        sep_r = reflct.split('_')
        if sep_h[4]==sep_r[4] and sep_h[5]==sep_r[5]:
            files_to_masaic_height.append(height)
            files_to_masaic_reflct.append(reflct)   

print(files_to_masaic_height, files_to_masaic_reflct)

#Defining a function to get extent of each tiff file
def get_extent(fn):
    "'returns min_x, max_y, max_x, min_y'"
    ds = gdal.Open(fn)
    gt = ds.GetGeoTransform()
    return (gt[0], gt[3], gt[0]+gt[1]*ds.RasterXSize, gt[3]+gt[5]*ds.RasterYSize)

#getting the extent of 1st file in the CHM file list 
min_x, max_y, max_x, min_y = get_extent(files_to_masaic_height[0])

#Getting the extent of mosaic file
for fn in files_to_masaic_height[1:]:
    minx, maxy, maxx, miny = get_extent(fn)
    min_x = min(min_x, minx)
    max_y = max(max_y, maxy)
    max_x = max(max_x, maxx)
    min_y = min(min_y, miny)

#Getting the number of rows and columns needed to creat mosaic file
in_ds = gdal.Open(files_to_masaic_height[0])
gt = in_ds.GetGeoTransform()
rows = math.ceil((max_y-min_y)/-gt[5])
columns = math.ceil((max_x-min_x)/gt[1])

#Creating a blank mosaic file
driver = gdal.GetDriverByName("GTiff")
out_ds = driver.Create('Height.tif', columns, rows)
out_ds.SetProjection(in_ds.GetProjection())
out_band = out_ds.GetRasterBand(1)

#georeferencing the mosaic file
gt = list (in_ds.GetGeoTransform())
gt[0], gt[3]=min_x, max_y
out_ds.SetGeoTransform(gt)
#reading each individual files and writing them insde the mosaic file
for fn in files_to_masaic_height:
    #getting the location where to write the data
    in_ds = gdal.Open(fn)
    gt = in_ds.GetGeoTransform()
    minx=gt[0]
    maxy=gt[3]
    x=math.ceil(minx-min_x)
    y=math.ceil(-(maxy-max_y))
    #reading input data as array
    data = in_ds.GetRasterBand(1).ReadAsArray()
    #writing the data 
    out_band.WriteArray(data,x,y)
del out_band, out_ds
"""at this point height data mosaic is complete. We will use the mosaiced height data extent to create another tiff file with 426 bands to mosaic hyperspectral image because our hyperspectral data has 426 bands"""
#creating 426 bands tiff file with same extent of height mosaic
chm= gdal.Open('Height.tif')
bn=chm.GetRasterBand(1)
driver = gdal.GetDriverByName("GTiff")
out = driver.Create('Reflectance.tif', bn.XSize, bn.YSize, 426, 2)
out.SetProjection(chm.GetProjection())
out.SetGeoTransform(chm.GetGeoTransform())

"""reading each reflectance file as numpy nd array. 426 reflectance bands create 426 dimentional numpy array. On the same location of that file inside the mosaic file there arrays are written according to their bands"""
for fn in files_to_masaic_reflct:
    #getting the location where to write the data
    hdf_ds = gdal.Open(fn)
    gt = hdf_ds.GetMetadata()['SRER_Reflectance_Reflectance_Data_Spatial_Extent_meters']
    sep =gt.split(' ')
    minx=float(sep[0])
    maxy=float(sep[3])
    x=math.ceil(minx-min_x)
    y=math.ceil(-(maxy-max_y))
    #creating input numpy nd array
    hdf_subds = hdf_ds.GetSubDatasets()
    reflct= hdf_subds[16][0]
    in_ds = gdal.Open(reflct)
    data = in_ds.ReadAsArray()
    #writing each array in each band
    for k,l in zip(list(range(1,427)),list(range(0,426))):
        out_band= out.GetRasterBand(k)
        in_array=data[:,:,l]
        out_band.WriteArray(in_array,x,y) 
del in_ds, out, out_band   





