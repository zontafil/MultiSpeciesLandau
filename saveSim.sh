#!/bin/sh

if [[ -z $1 ]]; then
    echo "Usage: saveSim.sh LABEL [PLOTDIR]"
    exit
fi

DATE=$(date +%Y%m%d_%H%M)
OUTFILE=$1_$DATE
cp main.cpp examples/$OUTFILE.cpp

if [[ -n $2 ]]; then
    cp out/out.avi $2/$OUTFILE.avi
fi