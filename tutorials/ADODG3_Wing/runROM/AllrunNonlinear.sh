#!/bin/bash

exec=mpirun
nProcs=4
solver=simpleROMFoam
runEndTime=1000
nSamples=90
nDVs=18
sampleCSV="shapeSamples.csv"
predictCSV="shapePredictions.csv"
predictSamples="1"

######################################################
# pre-processing
######################################################
 
# generate mesh
echo "Generating mesh.."
python genWingMesh.py > log.meshGeneration
plot3dToFoam -noBlank volumeMesh.xyz >> log.meshGeneration
autoPatch 60 -overwrite >> log.meshGeneration
createPatch -overwrite >> log.meshGeneration
renumberMesh -overwrite >> log.meshGeneration
echo "Generating mesh.. Done!"

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
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=train --nSamples=$nSamples --runEndTime=$runEndTime
  # run checkMesh for mesh quality
  $exec -np $nProcs checkMesh $pFlag > checkMeshLog
  # run the flow solver
  $exec -np $nProcs $solver $pFlag > flowLog
  cat objFuncs.dat
  cd ../runROM  

done

# deform the mesh to the reference point (last CFD sample)
$exec -np $nProcs python runFlow.py --task=deform --sample=$nSamples --mode=train --nSamples=$nSamples --runEndTime=$runEndTime
sleep 3

# link the simulations results from the CFD samples to runROM folder such that 
# the simpleROMFoam can read then to create the snapshot matrix
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
sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict

# run the offlineROM
$exec -np $nProcs $solver -mode offlineNonlinear $pFlag

# now we can clear the samples
for n in `seq 1 1 $nSamples`; do
    rm -rf processor*/$n
    rm -rf $n
done

######################################################
# runOnline
######################################################

# loop over all the prediction points
for n in $predictSamples; do

  # copy the runROM folder to prediction* folder
  rm -rf ../prediction$n
  cp -r ../runROM ../prediction$n
  cd ../prediction$n

  # copy the nearest sample field to 0 as initial field
  ((nM1=n-1))
  e=$(python findNearestSample.py $sampleCSV $predictCSV $nM1 2>&1)
  if [ $nProcs -gt 1 ]; then
    ((nProcsM1=nProcs-1))
    for p in `seq 0 1 $nProcsM1`; do
      rm -rf processor${p}/0/*
      cp -r ../sample$e/processor${p}/$runEndTime/* processor${p}/0/
    done
  else
    rm -rf 0/*
    cp -r ../sample$e/$runEndTime/* 0/
  fi

  # deform but not running the flow, now the geometry is predict sample but the based field is at initial field
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=predict --nSamples=$nSamples --runEndTime=$runEndTime

  # we want the onlineROM to use the latest CFD sample as the initial field
  sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict

  # run onlineROM, output UROM variables
  $exec -np $nProcs $solver -mode onlineNonlinear $pFlag

  # now run CFD at the predict sample
  echo "Run the flow at sample = $n"
  sed -i "/startFrom/c\startFrom       startTime;" system/controlDict
  sed -i "/solveAdjoint/c\solveAdjoint           false;" system/adjointDict
  $exec -np $nProcs $solver $pFlag > flowLog_${n}

  # print the result from CFD for reference
  cat objFuncs.dat

  cd ../runROM

done
