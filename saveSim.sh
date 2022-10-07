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
cp *.py plots/$OUTFILE/

cp out/*.png plots/$OUTFILE/png/
cp out/*.eps plots/$OUTFILE/eps/

mkdir -p savedData/$OUTFILE
cp out/data/*.txt savedData/$OUTFILE/

# save frames every 10th
cp $(printf '%s\n' out/data/*?????.png | awk 'NR%10 == 1') plots/$OUTFILE/frames/