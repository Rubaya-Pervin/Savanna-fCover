#!/usr/bin/env python
# coding: utf-8


import numpy as np
import glob, math, subprocess, h5py
from osgeo import gdal
import time
start = time.time()

#Read in reflectance hdf5 file 
hdf5_file = h5py.File('NEON_D14_SRER_DP3_511000_3520000_reflectance.h5','r')
file_attrs_string = str(list(hdf5_file.items()))
print(file_attrs_string)
file_attrs_string_split = file_attrs_string.split("'")
sitename = file_attrs_string_split[1]
refl = hdf5_file[sitename]['Reflectance']
#Extract the wavelength datasets
metadata = {}
metadata['wavelength'] = refl['Metadata']['Spectral_Data']['Wavelength'].value
#Extract bad band windows
metadata['bad_band_window1'] = (refl.attrs['Band_Window_1_Nanometers'])
metadata['bad_band_window2'] = (refl.attrs['Band_Window_2_Nanometers'])


data = gdal.Open('Reflct_red.tif')
data_clean = data.ReadAsArray()


#remove bad bands (NEON)
#1. define indices corresponding to min/max center wavelength for each bad band window:
bb1_ind0 = np.max(np.where((np.asarray(metadata['wavelength'])<float(metadata['bad_band_window1'][0]))))
bb1_ind1 = np.min(np.where((np.asarray(metadata['wavelength'])>float(metadata['bad_band_window1'][1]))))

bb2_ind0 = np.max(np.where((np.asarray(metadata['wavelength'])<float(metadata['bad_band_window2'][0]))))
bb2_ind1 = np.min(np.where((np.asarray(metadata['wavelength'])>float(metadata['bad_band_window2'][1]))))

bb3_ind0 = len(metadata['wavelength'])-10

#define valid band ranges from indices:
vb1 = list(range(0,bb1_ind0)); 
vb2 = list(range(bb1_ind1,bb2_ind0))
vb3 = list(range(bb2_ind1,bb3_ind0))

valid_band_range = [i for j in (range(0,bb1_ind0),
                                range(bb1_ind1,bb2_ind0),
                                range(bb2_ind1,bb3_ind0)) for i in j]

data_clean = data_clean[:,:,vb1+vb2+vb3]

metadata_clean['wavelength'] = [metadata['wavelength'][i] for i in valid_band_range]

#creating 426 bands tiff file with same extent of height mosaic
chm= gdal.Open('Height.tif')
bn=chm.GetRasterBand(1)
driver = gdal.GetDriverByName("GTiff")
out = driver.Create('Reflct_BigR.tif', bn.XSize, bn.YSize, 360, 2)
out.SetProjection(chm.GetProjection())
out.SetGeoTransform(chm.GetGeoTransform())
for k,l in zip(list(range(1,361)),list(range(0,360))):
	out.GetRasterBand(k).WriteArray(data[:,:,l]) 
del out
end = time.time()
print("tot time: ", (end-start))







