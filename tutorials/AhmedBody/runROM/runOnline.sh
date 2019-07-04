#!/bin/bash

exec=mpirun
nProcs=2

predictAngles="10 15 20"

predictI=1
for n in $predictAngles; do

  cp -r ../runROM ../prediction$predictI
  cd ../prediction$predictI

  $exec -np $nProcs python runFlow.py --task=write --angle=$n
  cp system/adjointDict.bk system/adjointDict
  cp system/controlDict.bk system/controlDict

  $exec -np $nProcs simpleFoamOnlineROM -parallel

  cd ../runROM

  ((predictI++))

done
