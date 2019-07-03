#!/bin/bash

p=1
angles="5 10 15 20 25"

# pre-processing
blockMesh
surfaceFeatureExtract
snappyHexMesh -overwrite
renumberMesh -overwrite
cp -r 0.orig 0


for n in $angles; do

  rm -rf ../sample$n
  cp -r ../runROM ../sample$n

  cd ../sample$n
  killall -9 foamRun.sh
  ./foamRun.sh $p &
  sleep 5
  mpirun -np $p python runFlow.py --angle $n
  killall -9 foamRun.sh

  cd ../runROM
  
done

python runFlow.py --task=write --angle=15.0

cp system/adjointDict.bk system/adjointDict
cp system/controlDict.bk system/controlDict

cp -r ../sample15/500 . 

offlineROM
onlineROM
