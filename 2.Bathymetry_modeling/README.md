Here are the primary scripts for constructing C-SHELPh commands:

Bathy\_utils.py: Contains all of the functions that C-SHELPh needs

run\_bathy.py: Calls the functions from Bathy\_utils.py in order to extract bathymetry photons

create\_run\_shell.py: A helper script to iterate over a range of threshold values and IS2 lasers. Creates run\_bathy.sh.
