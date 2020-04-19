.. _Installation:

Installation 
------------

ROMFoam depends on multiple prerequisites, to install ROMFoam

 - Install MACH framework following this website https://dafoam.readthedocs.io/en/latest/Installation.html **NOTE**: please follow the folder structure suggested in the installation documentation. Install all the dependencies except for DAFoam (no need to run its regression test either). This is because we have included DAFoam in ROMFoam. So DAFoam will be compiled along with ROMFoam later.
 
 - Download Slepc-3.6.3 at http://slepc.upv.es/download/distrib/slepc-3.6.3.tar.gz  Move the downloaded tarball to your $HOME/packages and untar it. Then go into the slepc-3.6.3 folder and run::

       ./configure && make SLEPC_DIR=$PWD PETSC_DIR=$HOME/packages/petsc-3.6.4 PETSC_ARCH=real-opt

 - Append the following to your $HOME/.bashrc file::
  
       # SLEPC
       export SLEPC_DIR=$HOME/packages/slepc-3.6.3
       export SLEPC_ARCH=real-opt
       export PATH=$SLEPC_DIR/$SLEPC_ARCH/bin:$PATH
       export PATH=$SLEPC_DIR/$SLEPC_ARCH/include:$PATH
       export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SLEPC_DIR/$SLEPC_ARCH/lib
       export SLEPC_LIB=$SLEPC_DIR/$SLEPC_ARCH/lib

   Then source it::

       source $HOME/.bashrc
 
 - Finally, download the ROMFoam repo, rename it to romfoam and put it in the $HOME/repos folder, load the OpenFOAM environment::
 
       source $HOME/OpenFOAM/OpenFOAM-v1812/etc/bashrc

   and go into the $HOME/repos/romfoam folder, and run::
   
       ./Allwmake


  