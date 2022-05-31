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

cp out/EnergySpecies.png plots/$OUTFILE/png/Energy.png
cp out/EnergyError.png plots/$OUTFILE/png/EnergySpecies.png
cp out/P.png plots/$OUTFILE/png/P.png
cp out/P.eps plots/$OUTFILE/eps/Vaxes.png
cp out/PerrNorm.png plots/$OUTFILE/png/PerrNorm.png
cp out/TemperatureSpecies.png plots/$OUTFILE/png/TemperatureSpecies.png
cp out/TemperatureSpeciesAxes.png plots/$OUTFILE/png/TemperatureSpeciesAxes.png

cp out/EnergySpecies.eps plots/$OUTFILE/eps/Energy.eps
cp out/EnergyError.eps plots/$OUTFILE/eps/EnergySpecies.eps
cp out/P.eps plots/$OUTFILE/eps/P.eps
cp out/P.eps plots/$OUTFILE/eps/Vaxes.eps
cp out/PerrNorm.eps plots/$OUTFILE/png/PerrNorm.eps
cp out/TemperatureSpecies.eps plots/$OUTFILE/eps/TemperatureSpecies.eps
cp out/TemperatureSpeciesAxes.eps plots/$OUTFILE/eps/TemperatureSpeciesAxes.eps

mkdir -p savedData/$OUTFILE
cp out/data/*.txt savedData/$OUTFILE/

cp $(printf '%s\n' out/data/*?????.png | awk 'NR%10 == 1') plots/$OUTFILE/frames/