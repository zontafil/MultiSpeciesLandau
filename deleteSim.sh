#!/bin/sh

if [ -z $1 ]; then
    echo "Usage: saveSim.sh examples/LABEL.cpp"
    exit
fi

SIMNAME=${1%.cpp}
SIMNAME="${SIMNAME##*/}"

rm examples/$SIMNAME.cpp
rm -rf plots/$SIMNAME
rm -rf savedData/$SIMNAME
