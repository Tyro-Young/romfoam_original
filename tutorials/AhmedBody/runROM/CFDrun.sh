#!/bin/bash

exec=mpirun
nProcs=4
solver=simpleROMFoam
runEndTime=500
avgFieldEvery=100
nSamples=5
resSamples=23
nDVs=1
predictSamples="1 2 3 4"

######################################################
# pre-processing
######################################################

# Generate the mesh, choose either way
#
##************** Option 1: Mesh generation using 1 CPU core ***********
blockMesh
surfaceFeatureExtract
snappyHexMesh -overwrite
renumberMesh -overwrite
./
##************** Option 2: Mesh generation using multiple CPU cores ***********
#blockMesh
#sed -i "/numberOfSubdomains/c\numberOfSubdomains $nProcs;" system/decomposeParDict
#surfaceFeatureExtract
#decomposePar
#$exec -np $nProcs snappyHexMesh -parallel
#reconstructParMesh -latestTime
#renumberMesh -overwrite
#rm -rf constant/polyMesh/*
#if [ -d "3/polyMesh" ]; then
#  mv 3/polyMesh/* constant/polyMesh/
#else
#  mv 2/polyMesh/* constant/polyMesh/
#fi
#rm -rf 1 2 3
#rm -rf processor*

# copy initial field and decompose the domain
cp -r 0.orig 0
sed -i "/numberOfSubdomains/c\numberOfSubdomains $nProcs;" system/decomposeParDict
if [ $nProcs -gt 1 ]; then
  decomposePar
fi

######################################################
# runOffline
######################################################

pFlag='-parallel'
if [ $nProcs -eq 1 ]; then
  pFlag=' '
fi

# Generate CFD samples
for n in `seq 1 1 $nSamples`; do

  rm -rf ../sample$n
  cp -r ../runROM ../sample$n

  cd ../sample$n
  # deform the mesh
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=train --nSamples=$nSamples --runEndTime=$runEndTime --avgFieldEvery=$avgFieldEvery
  # run checkMesh for mesh quality
  $exec -np $nProcs checkMesh $pFlag > checkMeshLog
  # run the flow solver
  $exec -np $nProcs $solver $pFlag > flowLog
  if [ $avgFieldEvery -gt 0 ]; then
    echo "Assigning mean to inst fields..."
    sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict
    $exec -np $nProcs meanToInstFields -varNames '(U p phi)'  $pFlag > meanToInstFieldsLog
    $exec -np $nProcs $solver -mode evalObj $pFlag > evalObjLog
  fi
  cat objFuncs.dat
  cd ../runROM

done

