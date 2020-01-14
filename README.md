ROMFoam: Reduced Order Modeling for OpenFOAM
===========================================

ROMFoam contains libraries and solvers to perform reduced-order modeling using OpenFOAM. 

Installation
------------

ROMFoam depends on multiple prerequisites, to install ROMFoam

 - Install MACH framework and all its dependencies following this website https://dafoam.readthedocs.io/en/latest/Installation.html **NOTE**: please follow the folder structure suggested in the installation documentation.
 
 - Download Slepc-3.6.3 at http://slepc.upv.es/download/distrib/slepc-3.6.3.tar.gz  Move the downloaded tarball to your $HOME/packages and untar it. Then go into the slepc-3.6.3 folder and run:

       ./configure && make SLEPC_DIR=$PWD PETSC_DIR=$HOME/packages/petsc-3.6.4 PETSC_ARCH=real-opt

 - Append the following to your $HOME/.bashrc file:
  
       # SLEPC
       export SLEPC_DIR=/home/ping/packages/slepc-3.6.3
       export SLEPC_ARCH=real-opt
       export PATH=$SLEPC_DIR/$SLEPC_ARCH/bin:$PATH
       export PATH=$SLEPC_DIR/$SLEPC_ARCH/include:$PATH
       export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SLEPC_DIR/$SLEPC_ARCH/lib
       export SLEPC_LIB=$SLEPC_DIR/$SLEPC_ARCH/lib
       # DAFoam dir
       export DAFOAM_DIR=$HOME/repos/dafoam/

   Then source it:

       source $HOME/.bashrc
 
 - Finally, download the OF_ROM repo, rename it to of_rom and put it in the $HOME/repos folder, load the OpenFOAM environment:
 
        source $HOME/OpenFOAM/OpenFOAM-v1812/etc/bashrc

   and go into the $HOME/repos/of_rom folder, and run::
   
       ./Allwmake

Tutorial
--------

Go to tutorial/AhmedBody/runROM, and Run:

    ./AllrunLinear.sh

The default setup has one design variable (ramp angle). It will first run the OpenFOAM's built in incompressible flow solver (simpleFoam) to generate samples with various ramp angle values (stored at /tutorial/AhmedBody/sample*). Then it will run simpleROMFoam to perform SVD and compute the reduced dR/dW and dR/dX matrices and save them to the disk. Finally, it will run simpleROMFoam to predict new fields given a new set of design variables.

To post-process the results, go to tutorial/AhmedBody/predict*. Open ParaView and load the paraview.foam file. The flow field variables computed by simpleFoam have names such as "U", "p", etc., while the flow fields computed from the reduced-order modeling have names such as "UROM", "pROM", etc. The modes from the SVD have names such as "UMode0", "UMode1", etc.

The setup have four major variations, and their design variables are: rampAngle, rideHeight, rampAngleAndRideHeight, and shape. The default is rampAngle, defined in --optVars in runFlow.py. You can select a different design variable setup, such as shape. In this case, you need to change --optVars to shape in runFlow.py, you also need to nSamples=20 and nDVs=4 in Allrun.sh. After this, run:

    ./AllrunLinear.sh

By default, each case does two prediction, the predicted point is defined by "DVs_Predict" in runFlow.py. If you want to do more prediction, change this parameter and also change predictSamples="1 2" in Allrun.sh. Here predictSamples="1 2" means prediction the first and second elements in "DVs_Predict".

For nonlinear ROM solvers, run:

    ./AllrunNonlinear.sh
