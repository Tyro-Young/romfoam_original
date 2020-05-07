#!/bin/bash

rm -rf 0
rm -rf postProcessing
rm -rf constant/extendedFeatureEdgeMesh
rm -rf constant/triSurface/*eMesh*
rm -rf *.bin *.info *Log* *.dat 
rm -rf jobFinished runCheckMesh* runFlowSolver* runAdjointSolver* runColoring
rm -rf constant/polyMesh
rm -rf processor*
rm -rf {1..9}*
rm -rf ../sample* ../prediction*
