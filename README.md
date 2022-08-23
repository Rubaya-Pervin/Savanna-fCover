# Savanna-fCover
This repository contains the data and scripts associated with Pervin et al. "Fusion of airborne hyperspectral and LiDAR canopy-height data for estimating fractional cover of tall woody plants, herbaceous vegetation, and other soil cover types in a semi-arid savanna ecosystem"

---------------------------------------------
## Files and Folders

code: code used in the analysis

  Merge.py : Merged data tiles for the whole study area to visually inspect presence of seam-lines and missing data.
  
  Data_Cleaning.py : Code to remove bad bands, no data values, and apply scale factor for hyperspectral data.
  
  BRDF correction is applied using codes from https://github.com/MarconiS/Estimating-individual-level-plant-traits-at-scale/
blob/master/py/hiperspectral_extractor2.py. Then, the following codes are applied to BRDF corrected data. 

  UnmixingBRDF_10T_sig.py : Codes to apply linear spectral unmixing on hyperspectral image only
  
  UnmixingBRDF_Height_10T_sig.py : fusion of hyperspectral and LiDAR CHM and linear unmixing of multimodal data. 
  
  fCover_Merge.py : merged fractional cover of different classes
