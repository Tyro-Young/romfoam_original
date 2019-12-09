#!/bin/bash

exec=mpirun
nProcs=1
solver=simpleROMFoam
runEndTime=1000
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

for n in `seq 1 1 $nSamples`; do

  rm -rf ../sample$n
  cp -r ../runROM ../sample$n

  cd ../sample$n
  killall -9 foamRun.sh
  ./foamRun.sh $exec $nProcs $solver &
  sleep 3
  $exec -np $nProcs python runFlow.py --task=run --sample=$n --mode=train --nSamples=$nSamples
  killall -9 foamRun.sh
  sleep 3

  cd ../runROM
  
done

$exec -np $nProcs python runFlow.py --task=writedelmat --sample=$nSamples --mode=train --nSamples=$nSamples
sleep 3
$exec -np $nProcs python runFlow.py --task=deform --sample=$nSamples --mode=train --nSamples=$nSamples
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

sed -i "/solveAdjoint/c\    solveAdjoint           true;" system/adjointDict
sed -i "/useColoring/c\    useColoring           true;" system/adjointDict
sed -i "/nFFDPoints/c\    nFFDPoints           $nDVs;" system/adjointDict
sed -i "/startFrom/c\    startFrom       latestTime;" system/controlDict

if [ $nProcs -eq 1 ]; then
  $solver -mode offlineLinear
else
  $exec -np $nProcs $solver -mode offlineLinear -parallel
fi

######################################################
# runOnline
######################################################

# calc refFields
rm -rf processor*
rm -rf {1..100}
killall -9 foamRun.sh
./foamRun.sh $exec $nProcs $solver &
sleep 3
$exec -np $nProcs python runFlow.py --task=run --sample=$refSample --mode=train --nSamples=$nSamples
killall -9 foamRun.sh
sleep 3


for n in $predictSamples; do

  rm -rf ../prediction$n
  cp -r ../runROM ../prediction$n
  cd ../prediction$n

  # deform but not running the flow, now the geometry is predict sample but the based field is at refSample
  $exec -np $nProcs python runFlow.py --task=deform --sample=$n --mode=predict --nSamples=$nSamples

  # run ROM, output UROM variables
  sed -i "/solveAdjoint/c\solveAdjoint           true;" system/adjointDict
  sed -i "/useColoring/c\useColoring           true;" system/adjointDict
  sed -i "/nFFDPoints/c\    nFFDPoints           $nDVs;" system/adjointDict
  sed -i "/startFrom/c\startFrom       latestTime;" system/controlDict
  if [ $nProcs -eq 1 ]; then
    $solver -mode onlineLinear
  else
    $exec -np $nProcs $solver -mode onlineLinear -parallel
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
