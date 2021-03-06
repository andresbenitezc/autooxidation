####
# Bash script to automate autooxidation Python scripts
# $File should be the file name with the spectrophotometer data
# $REFF is the folder containing the reference spectrums
# The number at the end references the column in the file for each replicate.
# Author: Andres Benitez
# Cite: Benitez Cardenas AS, Samuel PP, Olson JS. Current Challenges in the Development of Acellular Hemoglobin Oxygen Carriers by Protein Engineering. Shock. 2019;52(1S Suppl 1):28‐40. doi:10.1097/SHK.0000000000001053
####


#!/bin/sh

DIRECTORY="HbA 100uM"
SRC=$PWD

if [ ! -d "$DIRECTORY" ]
	 then 
		mkdir "$DIRECTORY"
fi
cd ./"$DIRECTORY"

REFF="$SRC/referencespectrums/ref_spectrum_HbApH5.csv"
FILE="$SRC/rawdata/170821_HbA100uMpH5.csv"
python3 ~/Tools/autoox-read.py $FILE $REFF 0 
python3 ~/Tools/autoox-read.py $FILE $REFF 1 
python3 ~/Tools/autoox-read.py $FILE $REFF 2 

REFF="$SRC/referencespectrums/ref_spectrum_HbApH6.csv"
FILE="$SRC/rawdata/170819_HbA100uMpH6.csv"
python3 ~/Tools/autoox-read.py $FILE $REFF 0 
python3 ~/Tools/autoox-read.py $FILE $REFF 1 
python3 ~/Tools/autoox-read.py $FILE $REFF 2 

REFF="$SRC/referencespectrums/ref_spectrum_HbA.csv"
FILE="$SRC/rawdata/170802_HbA100uM.csv"
python3 ~/Tools/autoox-read.py $FILE $REFF 0 
python3 ~/Tools/autoox-read.py $FILE $REFF 1 
python3 ~/Tools/autoox-read.py $FILE $REFF 2 

REFF="$SRC/referencespectrums/ref_spectrum_HbApH8.csv"
FILE="$SRC/rawdata/170807_HbA100uMpH8.csv"
python3 ~/Tools/autoox-read.py $FILE $REFF 0 
python3 ~/Tools/autoox-read.py $FILE $REFF 1 
python3 ~/Tools/autoox-read.py $FILE $REFF 2 

REFF="$SRC/referencespectrums/ref_spectrum_HbApH9.csv"
FILE="$SRC/rawdata/170811_HbA100uMpH9.csv"
python3 ~/Tools/autoox-read.py $FILE $REFF 0 
python3 ~/Tools/autoox-read.py $FILE $REFF 1 
python3 ~/Tools/autoox-read.py $FILE $REFF 2 

REFF="$SRC/referencespectrums/ref_spectrum_HbApH10.csv"
FILE="$SRC/rawdata/170816_HbA100uMpH10.csv"
python3 ~/Tools/autoox-read.py $FILE $REFF 0 
python3 ~/Tools/autoox-read.py $FILE $REFF 1 
python3 ~/Tools/autoox-read.py $FILE $REFF 2 
cd ../

