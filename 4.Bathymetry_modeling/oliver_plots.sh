#!/bin/bash

# change directories
cd ~/Desktop/summer_2021/icesat/data/raw/cuba

# create list of threshold values to loop through
declare -a threshold_array=("20" "25" "30" "35" "40" "45" "50" "55" "60" "65" "70" "75")

for filename in *.h5
do
	for value in $threshold_array
	do 
		python run.py -i $filename -l 3 -th $value
	done
done
