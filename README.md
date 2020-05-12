# autooxidation

R and Python scripts for processing spectrophotometer data for hemoglobin autooxidation experiments using spectral deconvolution to determine the relative concentrations of the different hemoglobin species in the sample over time.

There are 4 scripts included in the file for processing the data as well as an associated R package. 

* HbA100uM.sh is a bash script that helps automate multiple calls to the R and Python scripts.

* autoox-all-ggplot.R is an R script for plotting the output data using ggplot once the concentration of the different hemoglobin species have been determined.

* autoox-ratio.R is an R script for determing the hemoglobin species ratios over time after the data has been processed using the Python script autoox-read.py. The script uses the R package deconspectra to do the actual spectral decomposition.

* autoox-read.py reads the data from a .csv file, cleans up the data, and calls the R script autoox-ratio.R. The output of the R file is then written to an excel file and plotted.

* deconspectra is an R package that takes a mixed spectrum and any number of reference spectrums, and returns the ratios of each reference spectrum present in the mixed spectrum. 


R dependecies: scatterplot3d, minpack.lm, deconspectra (included in this repository), and signal.

Python dependecies: sys, re, numpy, and subprocess.
