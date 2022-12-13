# C-SHELPh

[![DOI](https://zenodo.org/badge/374786622.svg)](https://zenodo.org/badge/latestdoi/374786622)

A git repo to keep track of automatic retrieval of bathy photons and perform ML regression model

1.DownloadData: Provides an ICEPYX example (https://icepyx.readthedocs.io/en/latest/#) for automatically downloading IS2 data

2.Bathymetry\_modeling: Scripts required o run C-SHELPh

3.Regression\_Model: RSGISLib and Scikit-learn commands for performing Extra Trees Regression

4.Documents: Basic tutorial on using C-SHELPh

Note: The refraction correction currently uses photon level metrics (photon azimuth, photon elevation, satellite elev) to get an as accurate bathymetric model as possible. However, there is a small 'stepping' between IS2 segment id's. This is considered negligle given the quantities of data involved and value of the 'step' against the scale of the units
