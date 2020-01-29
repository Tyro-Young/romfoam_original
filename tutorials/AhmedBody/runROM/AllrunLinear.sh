#!/bin/bash

exec=mpirun
nProcs=2
solver=simpleROMFoam
runEndTime=500
nSamples=5
refSample=$nSamples
nDVs=1
predictSamples="1 2"

######################################################
# pre-processing
######################################################

# Generate the mesh
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

# Generate CFD samples
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

# write the dXv/dX matrix deltaVolMeshPionts.bin, which will be then used by the simpleROMFoam to compute dRdFFD
$exec -np $nProcs python runFlow.py --task=writedelmat --sample=$nSamples --mode=train --nSamples=$nSamples --runEndTime=$runEndTime
sleep 3

# deform the mesh to the reference point (last CFD sample)
$exec -np $nProcs python runFlow.py --task=deform --sample=$nSamples --mode=train --nSamples=$nSamples --runEndTime=$runEndTime
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
# modify the parameters in adjoiontDict and controlDict
sed -i "/solveAdjoint/c\    solveAdjoint           true;" system/adjointDict
sed -i "/useColoring/c\    useColoring           true;" system/adjointDict
sed -i "/nFFDPoints/c\    nFFDPoints           $nDVs;" system/adjointDict
sed -i "/startFrom/c\    startFrom       latestTime;" system/controlDict

# run the offlineROM
if [ $nProcs -eq 1 ]; then
  $solver -mode offlineLinear
else
  $exec -np $nProcs $solver -mode offlineLinear -parallel
fi

######################################################
# runOnline
######################################################

# clean up all the results, deform the mesh and run the flow at the refernece point
rm -rf processor*
for n in `seq 1 1 $nSamples`; do
    rm -rf $n
done
killall -9 foamRun.sh
./foamRun.sh $exec $nProcs $solver &
sleep 3
$exec -np $nProcs python runFlow.py --task=run --sample=$refSample --mode=train --nSamples=$nSamples --runEndTime=$runEndTime
killall -9 foamRun.sh
sleep 3

# loop over all the prediction points
for n in $predictSamples; do

  # copy the runROM folder to prediction* folder
  rm -rf ../prediction$n
  cp -r ../runROM ../prediction$n
  cd ../prediction$n

  # deform but not running the flow, now the geometry is predict sample but the based field is at refSample
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=predict --nSamples=$nSamples --runEndTime=$runEndTime

  # run onlineROM, output UROM variables
  sed -i "/solveAdjoint/c\solveAdjoint           true;" system/adjointDict
  sed -i "/useColoring/c\useColoring           true;" system/adjointDict
  sed -i "/nFFDPoints/c\    nFFDPoints           $nDVs;" system/adjointDict
  sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict
  if [ $nProcs -eq 1 ]; then
    $solver -mode onlineLinear
  else
    $exec -np $nProcs $solver -mode onlineLinear -parallel
  fi

  # now run CFD at the predict sample, overwrite the variable at refSample
  echo "Run the flow at sample = $n"
  sed -i "/startFrom/c\startFrom       startTime;" system/controlDict
  sed -i "/solveAdjoint/c\solveAdjoint           false;" system/adjointDict
  if [ $nProcs -eq 1 ]; then
    $solver > flowLog_${n}
  else
    $exec -np $nProcs $solver -parallel > flowLog_${n}
  fi
  
  # print the result from CFD for reference
  more objFuncs.dat

  killall -9 foamRun.sh

  cd ../runROM

done
