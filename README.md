# C-SHELPh V2

[![DOI](https://zenodo.org/badge/374786622.svg)](https://zenodo.org/badge/latestdoi/374786622)

A git repo to keep track of automatic retrieval of bathy photons and perform ML regression model

1.DownloadData: Provides an ICEPYX example (https://icepyx.readthedocs.io/en/latest/#) for automatically downloading IS2 data

2.Bathymetry\_modeling: Scripts required o run C-SHELPh

3.Regression\_Model: RSGISLib and Scikit-learn commands for performing Extra Trees Regression

4.Documents: Basic tutorial on using C-SHELPh

Note: The refraction correction currently uses photon level metrics (photon azimuth, photon elevation, satellite elev) to get an as accurate bathymetric model as possible. This requires some assuptions about the spacecraft orbit (its vertical movements are a constant rate thus linear interpolation holds true) and geometry and radius of the Earth. Ultimately, the advantages of increased precision of this implementation are believed to outweigh uncertainty introduced by the assumption. The refraction correction follows that outlined by the excellent Parrish et al., 2019 publication: https://www.mdpi.com/2072-4292/11/14/1634
