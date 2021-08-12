#!/bin/bash

# change directories
cd /Users/nmthomas/Documents/Research/ICESat2/SharkBay_ATL03

# create list of threshold values to loop through
declare -a threshold_array=(45)

for filename in *.h5
do
	for value in $threshold_array
	do 
		python /Users/nmthomas/Documents/Developer/Icesat2_bathymetry/4.Bathymetry_modeling/run.py -i $filename -l 1 -th $value
	done
done
