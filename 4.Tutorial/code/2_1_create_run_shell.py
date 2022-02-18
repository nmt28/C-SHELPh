import os
import subprocess
import glob

# Threshold values which control what density threshold is chosen to select the bathy photons. We have a selection of values here. For a full run we suggest using values 20-90 in interals of 5. This is something you will get a feel for as you use the code.
thresh_array = [5,10,15,20,25,30]

# The IS2 lasers. We will use 1 laser for now to reduce processing time but normally this would be [1,2,3]
laser_num = [1]

# The directory where the data is stored
dir = '../data/'

# grab all the files in the data directory where the file ends with .h5 (collect all the IS2 data)
files = glob.glob(dir + '/*.h5')

# create the file that you will run
txt_file = './2_2_run_bathy.sh'
# open the file and get ready to write to it
txt = open(txt_file, 'w+')

# iterate over the IS2 files
for file in files:
    # iterate over the lasers
    for l in laser_num:
        # iterate over the threshold values
        for t in thresh_array:
            # write the command to the shell file
            txt.write('python run_bathy.py -i ' + str(file) + ' -l ' + str(l) + ' -th ' + str(t) + '\n')
