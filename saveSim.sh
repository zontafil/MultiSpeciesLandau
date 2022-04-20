#!/bin/sh

if [ -z $1 ]; then
    echo "Usage: saveSim.sh LABEL"
    exit
fi

DATE=$(date +%Y%m%d_%H%M)
OUTFILE=${DATE}_$1

mkdir -p plots/$OUTFILE/eps
mkdir -p plots/$OUTFILE/png
mkdir -p plots/$OUTFILE/frames

cp main.cpp examples/$OUTFILE.cpp

cp out/out.avi plots/$OUTFILE/$OUTFILE.avi

cp out/EnergySpecies.png plots/$OUTFILE/png/${OUTFILE}_Energy.png
cp out/EnergyError.png plots/$OUTFILE/png/${OUTFILE}_EnergySpecies.png
cp out/minimumdist.png plots/$OUTFILE/png/${OUTFILE}_minimumdist.png
cp out/P.png plots/$OUTFILE/png/${OUTFILE}_P.png
cp out/PspeciesNorm.png plots/$OUTFILE/png/${OUTFILE}_PspeciesNorm.png
cp out/TemperatureSpecies.png plots/$OUTFILE/png/${OUTFILE}_TemperatureSpecies.png

cp out/EnergySpecies.eps plots/$OUTFILE/eps/${OUTFILE}_Energy.eps
cp out/EnergyError.eps plots/$OUTFILE/eps/${OUTFILE}_EnergySpecies.eps
cp out/minimumdist.eps plots/$OUTFILE/eps/${OUTFILE}_minimumdist.eps
cp out/P.eps plots/$OUTFILE/eps/${OUTFILE}_P.eps
cp out/PspeciesNorm.eps plots/$OUTFILE/eps/${OUTFILE}_PspeciesNorm.eps
cp out/TemperatureSpecies.eps plots/$OUTFILE/eps/${OUTFILE}_TemperatureSpecies.eps

mkdir -p savedData/$OUTFILE
cp out/data/*.txt savedData/$OUTFILE/

cp $(printf '%s\n' out/data/*?????.png | awk 'NR%10 == 1') plots/$OUTFILE/frames/