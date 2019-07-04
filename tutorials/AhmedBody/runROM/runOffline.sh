#!/bin/bash

exec=mpirun
nProcs=2
solver=simpleDAFoam
runEndTime=500
refAngle=15
angles="10 20 15"
nSamples=3

# pre-processing
blockMesh
surfaceFeatureExtract
snappyHexMesh -overwrite
renumberMesh -overwrite
cp -r 0.orig 0
sed -i "/numberOfSubdomains/c\numberOfSubdomains $nProcs;" system/decomposeParDict
decomposePar

sampleI=1
for n in $angles; do

  rm -rf ../sample$sampleI
  cp -r ../runROM ../sample$sampleI

  cd ../sample$sampleI
  killall -9 foamRun.sh
  ./foamRun.sh $exec $nProcs $solver &
  sleep 5
  $exec -np $nProcs python runFlow.py --angle=$n
  killall -9 foamRun.sh

  ((sampleI++))

  cd ../runROM
  
done

$exec -np $nProcs python runFlow.py --task=write --angle=$refAngle

cp system/adjointDict.bk system/adjointDict
cp system/controlDict.bk system/controlDict

((nProcsM1=nProcs-1))
for m in `seq 1 1 $nSamples`; do
  for n in `seq 0 1 $nProcsM1`; do
    cd processor${n}
    ln -s ../../sample${m}/processor${n}/$runEndTime $m
    cd ../
  done
done

$exec -np $nProcs simpleFoamOfflineROM -parallel

