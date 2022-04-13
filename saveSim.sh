#!/bin/sh

if [ -z $1 ]; then
    echo "Usage: saveSim.sh LABEL"
    exit
fi

DATE=$(date +%Y%m%d_%H%M)
OUTFILE=$1_$DATE
cp main.cpp examples/$OUTFILE.cpp

cp out/out.avi plots/$OUTFILE.avi
cp out/EnergySpecies.png plots/${OUTFILE}_Energy.png
cp out/EnergyError.png plots/${OUTFILE}_EnergySpecies.png