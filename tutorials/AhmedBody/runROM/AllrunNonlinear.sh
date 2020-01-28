#!/bin/bash

exec=mpirun
nProcs=1
solver=simpleROMFoam
runEndTime=500
nSamples=5
refSample=$nSamples
nDVs=1
predictSamples="1 2"

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

for n in `seq 1 1 $nSamples`; do

  rm -rf ../sample$n
  cp -r ../runROM ../sample$n

  cd ../sample$n
  killall -9 foamRun.sh
  ./foamRun.sh $exec $nProcs $solver &
  sleep 3
  $exec -np $nProcs python runFlow.py --task=run --sample=$n --mode=train --nSamples=$nSamples --runEndTime=$runEndTime
  killall -9 foamRun.sh
  sleep 3

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
if [ $nProcs -eq 1 ]; then
  $solver -mode offlineNonlinear
else
  $exec -np $nProcs $solver -mode offlineNonlinear -parallel
fi

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
rm -rf processor*/{1..100}
rm -rf {1..100}
killall -9 foamRun.sh

for n in $predictSamples; do

  rm -rf ../prediction$n
  cp -r ../runROM ../prediction$n
  cd ../prediction$n

  # deform but not running the flow, now the geometry is predict sample but the based field is at refSample
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=predict --nSamples=$nSamples --runEndTime=$runEndTime

  # run ROM, output UROM variables
  sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict
  if [ $nProcs -eq 1 ]; then
    $solver -mode onlineNonlinear
  else
    $exec -np $nProcs $solver -mode onlineNonlinear -parallel
  fi


  # now run the flow at the predict sample, overwrite the variable at refSample
  echo "Run the flow at sample = $n"
  sed -i "/startFrom/c\startFrom       startTime;" system/controlDict
  sed -i "/solveAdjoint/c\solveAdjoint           false;" system/adjointDict
  if [ $nProcs -eq 1 ]; then
    $solver > flowLog_${n}
  else
    $exec -np $nProcs $solver -parallel > flowLog_${n}
  fi
  more objFuncs.dat

  cd ../runROM

done
