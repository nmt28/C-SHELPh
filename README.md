# C-SHELPh V2

[![DOI](https://zenodo.org/badge/374786622.svg)](https://zenodo.org/badge/latestdoi/374786622)

### C-SHELPh is the Classification of Sub-aquatic Height Extracted Photons. It is designed to isolate bathymetric photons in ICESat2 ATL03 files.

## Installation

### It is recommended that the dependancies are installed via:
```
conda install -c conda-forge geopandas utm numpy matplotlib s3fs xarray zarr pyproj proj-data h5py

pip install cshelph

```

## Using C-SHELPh

### This software is accompanied by two run scripts:

* run_bathymetry_extraction.ipynb: A notebook designed to process single tracks to gain familiarity with the software
* run_bathy_extraction.py: Automated runs of c-shelph for mass processing

### A simple use of C-SHELPh is:
```
python run_bathy_extraction.py -i icesat2_atl03_file.h5 -l 1 -th 20
```
### where:
```
-i: the input ICESat2 ATL03 h5 file
-l: laser number (1-3)
-th: density threshold value (percentile; 0-100) which is used to change the sensitivity of the photon classification to noise
```
### Additional options can be specifed to customize runs and override defaults (which are based on some underlying assumptions) 

## FAQs

### **Why does C-SHELPh classify my water surface as bathy?**

### C-SHELPh performs an orthometric correction that models the water surface at approximately 0m elevation in reference to the EGM2008 geoid. In some locations this is not true and so C-SHELPh gets confused and thresholds set in the code to help C-SHELPh find the water surface are not sufficient. These can be changed in the code.

### **Why isn't C-SHELPh detecting my bathy points**

### C-SHELPh is classifies photons based on photon density. In cases where the noise is higher than the signal (e.g., a small quantity of steep/deep and/or turbid waters) C-SHELph struggles to perform well. Optimum conditions for C-SHELPh are clear shallow waters.

## Additional notes

### Note: The refraction correction currently uses photon level metrics (photon azimuth, photon elevation, satellite elev) to get an as accurate bathymetric model as possible. This requires some assuptions about the spacecraft orbit (its vertical movements are a constant rate thus linear interpolation holds true) and geometry and radius of the Earth. Ultimately, the advantages of increased precision of this implementation are believed to outweigh uncertainty introduced by the assumption. The refraction correction follows that outlined by the excellent Parrish et al., 2019 publication: https://www.mdpi.com/2072-4292/11/14/1634
