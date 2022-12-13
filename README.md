# C-SHELPh

[![DOI](https://zenodo.org/badge/374786622.svg)](https://zenodo.org/badge/latestdoi/374786622)

A git repo to keep track of automatic retrieval of bathy photons and perform ML regression model

1.DownloadData: Provides an ICEPYX example (https://icepyx.readthedocs.io/en/latest/#) for automatically downloading IS2 data

2.Bathymetry\_modeling: Scripts required o run C-SHELPh

3.Regression\_Model: RSGISLib and Scikit-learn commands for performing Extra Trees Regression

4.Documents: Basic tutorial on using C-SHELPh

Note: The refraction correction currently uses a simplification of the incidence angle (theta1) as outlined by Parrish et al, 2019. For greater precision, the curvature of the Earth needs to be accounted for along the ICESat-2 track. This will be integrated into C-SHELPh as time permits.
