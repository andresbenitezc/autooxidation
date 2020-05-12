####
# Python script to clean data and run R script for spectral deconvolution.
# The script requires three inputs; a file containing the spectrum to be analyzed, a file containing the reference spectrums,
# and a the number indicating which spectrum to analyze (the code was written to analyze data from triplicate experiments)
# Author: Andres Benitez
# Bissé E, Schaeffer-Reiss C, Van Dorsselaer A, et al. Hemoglobin Kirklareli (α H58L), a New Variant Associated with Iron Deficiency and Increased CO Binding. J Biol Chem. 2017;292(6):2542‐2555. doi:10.1074/jbc.M116.764274
####

#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

import sys
import re
import numpy as np
import subprocess

# Define a function to find the closest value in an array to a given value.
def closest(array, value):
    index = (np.abs(array - value)).argmin()
    return index

# The script arguments are the filename, which is the csv output from the spec, the sample number, and the name of the reference file containing the 
#reference spectrums for the given sample.
filename = sys.argv[1]
fileref = sys.argv[2]
sample = int(sys.argv[3])

# The time steps will be taken from a prepared file with the same name as the reference sample, in a folder called "time"
filetime = filename.replace('rawdata', 'time').replace('.csv', '-time.csv')


# Take specdata file and convert to a list of lists where each row is a list. Do the same with the time file
specdata = [[item.strip() for item in line.split(',')] for line in open(filename)]
timedata = [[float(item.strip()) for item in line.split(',')] for line in open(filetime)]
refdata = [[item.strip() for item in line.split(',')] for line in open(fileref)]

# Create time range
time = np.concatenate((np.arange(0, timedata[1][0], timedata[0][0]), 
	np.arange(timedata[1][0], timedata[1][1], timedata[0][1]), 
	np.arange(timedata[1][1], timedata[1][2] + timedata[0][2], timedata[0][2]))).tolist()

# Remove empty last row. 
specdata.pop()

# Use the names of the specdata columns to determine how many samples are in the file, and what their names are.
names = [name.replace('_1','') for name in specdata[0][0::2] if re.search('_1$', name)]

# Covert the list of lists from row to column (transverse)
specdata = list(map(list, zip(*specdata)))

# Remove last row, which is empty, and remove the first two rows, since they are not useful after removing the names
specdata.pop()
[item.pop(0) for item in specdata]
[item.pop(0) for item in specdata]

# The first list (column) is the wavelength column - all wavelength columns are equal
wavelength = specdata[0]

# Remove wavelength columns from the rest of the array
specdata = specdata[1::2]

# Subset for current sample being analyzed
specdata = specdata[sample::len(names)]

# Convert reference arrays to match wavelengths of experiment data, and extract names
refarray = np.asarray(list(map(float, refdata[0][1:])))
refindex = [closest(refarray, float(wave)) for wave in wavelength]

refnames = [row[0] for row in refdata][1:]
[item.pop(0) for item in refdata]

reffinal = [ [refdata[i][j] for j in refindex] for i in range(1,len(refdata)) ]

# Remove columns beyond the time that is in the time input file
specdata = specdata[0:len(time)]

# Insert wavelenth rows and reference rows
specdata.insert(0, wavelength)
[specdata.insert(1, row) for row in reversed(reffinal)]

# Convert from columns back to rows and write to file, and add a header including the time.
specrow = list(zip(*specdata))

outputfilename = re.sub('.+/rawdata/', '', filename.replace('.csv', '_') + names[sample] + '.csv')

outputfile = open(outputfilename, 'w+')
specrow.insert(0, ['Wavelength'] + refnames + [str(item) for item in time])

outputfile.write('\n'.join([','.join(row) for row in specrow]))
outputfile.close()

command = 'Rscript'
path2script = '/Users/olsonlab/Tools/autoox-ratio.R'
args = [outputfilename]

cmd = [command, path2script] + args

try:
	x = subprocess.check_output(cmd, universal_newlines = True)
except subprocess.CalledProcessError as e:
	x = e.output


