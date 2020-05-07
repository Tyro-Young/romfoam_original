#!/bin/bash

exec=mpirun
nProcs=1
solver=simpleROMFoam
runEndTime=1000
avgFieldEvery=200
nSamples=10
refSample=$nSamples
nDVs=2
predictSamples="1 2"

######################################################
# pre-processing
######################################################

cp -r 0.orig 0
sed -i "/numberOfSubdomains/c\numberOfSubdomains $nProcs;" system/decomposeParDict
if [ $nProcs -gt 1 ]; then
  decomposePar
fi

#####################################################
# runOffline
#####################################################

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

# deform the mesh to the reference point (last CFD sample)
$exec -np $nProcs python runFlow.py --task=deform --sample=$nSamples --mode=train --nSamples=$nSamples --runEndTime=$runEndTime --avgFieldEvery=$avgFieldEvery
sleep 3

# link the simulations results from the CFD samples to runROM folder such that 
# the simpleROMFoam can read then  to create the snapshot matrix
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

# modify the parameters in and controlDict because we want to use the 
# flow field from the latest sample to compute the preconditioner matrix
sed -i "/startFrom/c\    startFrom       latestTime;" system/controlDict

# run the offlineROM
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

# loop over all the prediction points
for n in $predictSamples; do

  # copy the runROM folder to prediction* folder
  rm -rf ../prediction$n
  cp -r ../runROM ../prediction$n
  cd ../prediction$n

  # deform but not running the flow, now the geometry is predict sample but the based field is at refSample
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=predict --nSamples=$nSamples --runEndTime=$runEndTime --avgFieldEvery=$avgFieldEvery

  # we want the onlineROM to use the latest CFD sample as the initial field
  sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict

  # run onlineROM, output UROM variables
  $exec -np $nProcs $solver -mode onlineNonlinear $pFlag

  # now run CFD at the predict sample
  echo "Run the flow at sample = $n"
  sed -i "/startFrom/c\startFrom       startTime;" system/controlDict
  sed -i "/solveAdjoint/c\solveAdjoint           false;" system/adjointDict
  $exec -np $nProcs $solver $pFlag > flowLog_${n}
  if [ $avgFieldEvery -gt 0 ]; then
    echo "Assigning mean to inst fields..."
    sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict
    $exec -np $nProcs meanToInstFields -varNames '(U p phi)'  $pFlag > meanToInstFieldsLog
    $exec -np $nProcs $solver -mode evalObj $pFlag > evalObjLog
  fi

  # print the result from CFD for reference
  cat objFuncs.dat

  cd ../runROM

done
