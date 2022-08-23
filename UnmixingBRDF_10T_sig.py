#!/N/u/macblab/Carbonate/.conda/envs/py37_env/bin/python
# coding: utf-8

import gdal
import h5py, random, glob, math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg') 
import numpy as np
import pandas as pd
#for fractional shrub coverage
import pysptools.util as util
import pysptools.eea as eea #endmembers extraction algorithms
import pysptools.abundance_maps as amap
import pysptools.classification as cls
import pysptools.material_count as cnt
import time
start = time.time()

def get_wavelength(h5,wd):
	hdf5_file = h5py.File(h5, 'r')
	file_attrs_string = str(list(hdf5_file.items()))
	file_attrs_string_split = file_attrs_string.split("'")
	sitename = file_attrs_string_split[1]
	refl = hdf5_file[sitename]['Reflectance']
	#Create dictionary containing relevant metadata information(from NEON algorithm)
	metadata = {}
	metadata['map info'] = refl['Metadata']['Coordinate_System']['Map_Info'].value
	metadata['wavelength'] = refl['Metadata']['Spectral_Data']['Wavelength'].value
	#Extract bad band windows
	metadata['bad_band_window1'] = (refl.attrs['Band_Window_1_Nanometers'])
	metadata['bad_band_window2'] = (refl.attrs['Band_Window_2_Nanometers'])
	metadata_clean= metadata.copy()
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
	wavelength = [metadata['wavelength'][i] for i in valid_band_range]
	return(wavelength)
	
def get_data(brdf_rfl):
	dt = gdal.Open(brdf_rfl)
	sep = brdf_rfl.split('_')
	refl = dt.ReadAsArray()
	data = np.transpose(refl, (1,2,0))
	return(data)

def get_sig4(data,wavelength,wd):
#	dt = gdal.Open(brdf_rfl)
#	sep_r = brdf_rfl.split('_')
#	refl = dt.ReadAsArray()
#	data = np.transpose(refl, (1,2,0))

	#endmember extraction algorithm
	ee = eea.NFINDR()
	U = ee.extract(data,4,maxit=4, normalize=False,ATGP_init=True)
	path = wd
	name = path + "10Tile_sig4."
	np.save(name + 'npy', U)

	plt.figure()
	plt.plot(wavelength, U[0,:], label='EM1')  
	plt.plot(wavelength, U[1,:], label='EM2')  
	plt.plot(wavelength, U[2,:], label='EM3')
	plt.plot(wavelength, U[3,:], label='EM4')

	plt.xlabel('Wavelength nm+height')
	plt.ylabel('Reflectance')
	plt.title("4 endmembers Spectral Signature with height: 10 tile")
	xticks = np.arange(0, 2501, 500)
	plt.xticks(xticks,('0','500', '1000', '1500', '2000', 'height'))
	plt.legend()
	plt.savefig(name + 'png')
	plt.show()
	return(U)

def get_sig5(data,wavelength,wd):
#	dt = gdal.Open(brdf_rfl)
#	sep_r = brdf_rfl.split('_')
#	refl = dt.ReadAsArray()
#	data = np.transpose(refl, (1,2,0))

	#endmember extraction algorithm
	ee = eea.NFINDR()
	U = ee.extract(data,5,maxit=5, normalize=False,ATGP_init=True)
	path = wd
	name = path + "10Tile_sig5."
	np.save(name + 'npy', U)

	plt.figure()
	plt.plot(wavelength, U[0,:], label='EM1')  
	plt.plot(wavelength, U[1,:], label='EM2')  
	plt.plot(wavelength, U[2,:], label='EM3')
	plt.plot(wavelength, U[3,:], label='EM4')
	plt.plot(wavelength, U[4,:], label='EM5')
	plt.xlabel('Wavelength nm+height')
	plt.ylabel('Reflectance')
	plt.title("5 endmembers Spectral Signature with height: 10 tile")
	xticks = np.arange(0, 2501, 500)
	plt.xticks(xticks,('0','500', '1000', '1500', '2000', 'height'))
	plt.legend()
	plt.savefig(name + 'png')
	plt.show()
	return(U)



def get_FC(brdf_rfl,U):
	dt = gdal.Open(brdf_rfl)
	refl = dt.ReadAsArray()
	data = np.transpose(refl, (1,2,0))
	#define am object using the amap
	am = amap.FCLS() 
	#create abundance maps for the HSI cubems
	amaps = am.map(data,U,normalize=False) 
	#diaplay
	#am.display(colorMap='jet',columns=4,suffix='SRER')
	return(amaps)

def save_FC(brdf_rfl,FC, wd):
	sep = brdf_rfl.split('_')
	hfile = '/N/project/SRER/CanopyHeightModelGtifReprocessed/NEON_D14_SRER_DP3_' +sep[1]+'_'+sep[2].split('.')[0]+ '_CHM.tif'
	hm = gdal.Open(hfile)
	b=hm.GetRasterBand(1)
	for k in list(range(0,FC.shape[2])):
		path = wd + "FC/" 
		file_name = path +'SRER_'+sep[1]+'_'+sep[2].split('.')[0]+'_cl'+np.str(FC.shape[2])+'_EM'+str(k+1)+'.tif'
		drive = gdal.GetDriverByName("GTiff")
		out = drive.Create(file_name, b.XSize, b.YSize, 1, b.DataType)
		out.SetProjection(hm.GetProjection())
		out.SetGeoTransform(hm.GetGeoTransform())
		out_band= out.GetRasterBand(1)
		in_array=FC[:,:,k]
		out_band.WriteArray(in_array)
		del out, out_band
		

#Getting list of  files 
h5_files = glob.glob('/N/project/SRER/Data/*.h5')
files_to_classify = glob.glob('/N/project/SRER/BRDFCorrected/*.tif')
files_to_get_sig = []
for reflct in files_to_classify:
	sep_r = reflct.split('_')
	if int(sep_r[1])!= 511000:
		files_to_get_sig.append(reflct)
random.seed(55)
sig_files = random.sample(files_to_get_sig,10)

for x in sig_files:
	print(x)
h5 = random.choice(h5_files)
print(h5)

#reading 10 tiles and combining them 
data = get_data(sig_files[0])
for x in list(range(1,10)):
	dt = get_data(sig_files[x])
	data = np.concatenate((data,dt), axis = 0)
print(data.shape)
wd = '/N/project/SRER/Unmixing/BRDF/10TSig/'
wavelength = get_wavelength(h5,wd)
#U4 = get_sig4(data,wavelength,wd)
U5 = get_sig5(data,wavelength,wd)
#classifying all tiles
for refl in files_to_classify:
	FC = get_FC(refl,U5)
	save_FC(refl,FC, wd)
end = time.time()
print("tot time: ", (end-start))
