#!/bin/bash

# pre-processing
blockMesh
surfaceFeatureExtract
snappyHexMesh -overwrite
renumberMesh -overwrite
cp -r 0.orig 0

offlineROM

# these are the actually commands to run the case
#./foamRun.sh &
#mpirun -np 2 python runScript.py
