#!/bin/bash

exec=mpirun
nProcs=2
solver=simpleROMFoam
runEndTime=500
nSamples=5
refSample=$nSamples
nDVs=1
predictSamples="1"

######################################################
# pre-processing
######################################################

blockMesh
surfaceFeatureExtract
snappyHexMesh -overwrite
renumberMesh -overwrite
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
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=train --nSamples=$nSamples --runEndTime=$runEndTime
  # run checkMesh for mesh quality
  $exec -np $nProcs checkMesh $pFlag > checkMeshLog
  # run the flow solver
  $exec -np $nProcs $solver $pFlag > flowLog
  cat objFuncs.dat
  cd ../runROM
  
done

$exec -np $nProcs python runFlow.py --task=deform --sample=$nSamples --mode=train --nSamples=$nSamples --runEndTime=$runEndTime
sleep 3

if [ $nProcs -gt 1 ]; then
  ((nProcsM1=nProcs-1))
  for m in `seq 1 1 $nSamples`; do
    for n in `seq 0 1 $nProcsM1`; do
      cd processor${n}
      ln -s ../../sample${m}/processor${n}/$runEndTime $m
      cd ../
    done
  done
else
  for m in `seq 1 1 $nSamples`; do
    ln -s ../sample${m}/$runEndTime $m
  done
fi

sed -i "/startFrom/c\    startFrom       latestTime;" system/controlDict

# use the last sample field as ref
$exec -np $nProcs $solver -mode offlineNonlinear $pFlag

######################################################
# runOnline
######################################################

# copy the last sample field to 0
if [ $nProcs -gt 1 ]; then
  ((nProcsM1=nProcs-1))
  for n in `seq 0 1 $nProcsM1`; do
    rm -rf processor${n}/0/*
    cp -r processor${n}/$nSamples/* processor${n}/0/
  done
else
  rm -rf 0/*
  cp -r $nSamples/* 0/
fi

# now we can clear the samples
for n in `seq 1 1 $nSamples`; do
    rm -rf processor*/$n
    rm -rf $n
done

for n in $predictSamples; do

  rm -rf ../prediction$n
  cp -r ../runROM ../prediction$n
  cd ../prediction$n

  # deform but not running the flow, now the geometry is predict sample but the based field is at refSample
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=predict --nSamples=$nSamples --runEndTime=$runEndTime

  # run ROM, output UROM variables
  sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict
  $exec -np $nProcs $solver -mode onlineNonlinear $pFlag
  cp objFuncs.dat ../runROM
  #echo "CD: 0.402261736832066 (ROM Ref)"


  # now run the flow at the predict sample, overwrite the variable at refSample
  #echo "Run the flow at sample = $n"
  #sed -i "/startFrom/c\startFrom       startTime;" system/controlDict
  #sed -i "/solveAdjoint/c\solveAdjoint           false;" system/adjointDict
  #$exec -np $nProcs $solver $pFlag > flowLog_${n}
  #more objFuncs.dat
  #echo "CD 0.4029108949571894 (ROM Ref)"

  cd ../runROM

done
