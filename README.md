OF_ROM: Reduced Order Modeling for OpenFOAM
===========================================

OF_ROM contains libraries and solvers to perform reduced-order modeling using OpenFOAM. 

Installation
------------

OF_ROM depends on multiple prerequisites, to install OF_ROM

 - Install DAFoam and all its dependences following this website https://dafoam.readthedocs.io/en/latest/Installation.html **NOTE**: please follow the folder structure suggested in the installation documentation.
 
 - Download Slepc-3.6.3 at http://slepc.upv.es/download/distrib/slepc-3.6.3.tar.gz  Move the downloaded tarball to your $HOME/packages and untar it. Then go into the slepc-3.6.3 folder and run::

    ./configure && make SLEPC_DIR=$PWD PETSC_DIR=$HOME/packages/petsc-3.6.4 PETSC_ARCH=real-opt

 - Append the following to your $HOME/.bashrc file::
  
    # SLEPC
    export SLEPC_DIR=/home/ping/packages/slepc-3.6.3
    export SLEPC_ARCH=real-opt
    export PATH=$SLEPC_DIR/$SLEPC_ARCH/bin:$PATH
    export PATH=$SLEPC_DIR/$SLEPC_ARCH/include:$PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SLEPC_DIR/$SLEPC_ARCH/lib
    export SLEPC_LIB=$SLEPC_DIR/$SLEPC_ARCH/lib
    # DAFoam dir
    export DAFOAM_DIR=$HOME/repos/dafoam/

   Then source it::

    source $HOME/.bashrc
 
 - Finally, download the OF_ROM repo, rename it to of_rom and put it in the $HOME/repos folder, load the OpenFOAM environment::
 
     source $HOME/OpenFOAM/OpenFOAM-v1812/etc/bashrc

   and go into the $HOME/repos/of_rom folder, and run::
   
    ./Allwmake

Tutorial
--------

Run::

 ./Allrun.sh

