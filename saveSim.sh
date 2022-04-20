#!/bin/sh

if [ -z $1 ]; then
    echo "Usage: saveSim.sh LABEL"
    exit
fi

DATE=$(date +%Y%m%d_%H%M)
OUTFILE=${DATE}_$1

mkdir -p plots/eps
mkdir -p plots/png

cp main.cpp examples/$OUTFILE.cpp

cp out/out.avi plots/$OUTFILE.avi

cp out/EnergySpecies.png plots/png/${OUTFILE}_Energy.png
cp out/EnergyError.png plots/png/${OUTFILE}_EnergySpecies.png
cp out/minimumdist.png plots/png/${OUTFILE}_minimumdist.png
cp out/P.png plots/png/${OUTFILE}_P.png
cp out/PspeciesNorm.png plots/png/${OUTFILE}_PspeciesNorm.png
cp out/TemperatureSpecies.png plots/png/${OUTFILE}_TemperatureSpecies.png

cp out/EnergySpecies.eps plots/eps/${OUTFILE}_Energy.eps
cp out/EnergyError.eps plots/eps/${OUTFILE}_EnergySpecies.eps
cp out/minimumdist.eps plots/eps/${OUTFILE}_minimumdist.eps
cp out/P.eps plots/eps/${OUTFILE}_P.eps
cp out/PspeciesNorm.eps plots/eps/${OUTFILE}_PspeciesNorm.eps
cp out/TemperatureSpecies.eps plots/eps/${OUTFILE}_TemperatureSpecies.eps

mkdir -p savedData/$OUTFILE
cp out/data/*.txt savedData/$OUTFILE/