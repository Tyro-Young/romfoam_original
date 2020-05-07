.. _Installation:

Installation 
------------

There are two options to install ROMFoam: **pre-compiled package** and **source code**. If you are running ROMFoam for the first time, we recommend using the pre-compiled version, which supports Linux (Ubuntu, Fedora, CentOS, etc), MacOS, and Windows systems. For production runs on an HPC system, you need to compile ROMFoam from the source.

#. **Pre-compiled package**

   The pre-compiled package is available on Docker Hub. Before downloading the pre-compiled package, you need to install **Docker**. Follow the installation instructions for `Ubuntu <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_, `Fedora <https://docs.docker.com/install/linux/docker-ce/fedora/>`_, `CentOS <https://docs.docker.com/install/linux/docker-ce/centos/>`_, `MacOS <https://docs.docker.com/docker-for-mac/install/>`_, and  `Windows <https://docs.docker.com/docker-for-windows/install/>`_. 
 
   For example, on Ubuntu 18.04, you can install the latest Docker by running this command in the terminal::

    sudo apt-get remove docker docker-engine docker.io containerd runc && sudo apt-get update && sudo apt-get install apt-transport-https ca-certificates curl gnupg-agent software-properties-common -y && curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add - && sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" && sudo apt-get update && sudo apt-get install docker-ce -y

   Then you need to add your user name to the Docker group by::

    sudo usermod -aG docker $USER

   After this, you need to **logout and re-login** your account to make the usermod command effective. Once done, verify the Docker installation by::

    docker --version

   You should be able to see your installed Docker version. Note that different operating systems have very different Docker installation process, refer to the above links for more details. 

   Once the Docker is installed and verified, run this command from the terminal::

    docker run -it --rm -u dafoamuser -v $HOME:/home/dafoamuser/mount -w /home/dafoamuser/mount dafoam/opt-packages:latest bash

   It will first download the pre-compiled package from the Docker Hub if it has not been downloaded. Then it will start a Docker container (a light-weight virtual machine), mount your local computer's home directory to the container's ``mount`` directory, login to ``mount`` as dafoamuser, and set the relevant environmental variables. 

   You should be in a Docker container environment now. Next, go to the repos folder and download ROMFoam::

    cd $HOME/repos && \
    git clone https://github.com/mdolab/romfoam && \
    cd romfoam && \
    ./Allwmake
   
   Now, you are ready to run a case. We recommend you copy the tutorial from ``$HOME/repos/romfoam/tutorials`` to ``$HOME/mount/tutorials`` and run jobs from ``$HOME/mount/tutorials``. This is because if you directly run jobs on $HOME/repos/romfoam/tutorials, the results **WILL BE DELETED** when you exit the Docker container. 

   When you are done, exit the Docker container by running::
   
    exit

   **NOTE:** Because we use Docker, everytime you want to run a ROMFoam job, you need to login to the Docker container, download the latest ROMFoam repo and compile it, i.e., doing this::

    docker run -it --rm -u dafoamuser -v $HOME:/home/dafoamuser/mount -w /home/dafoamuser/mount dafoam/opt-packages:latest bash && \
    cd $HOME/repos && \
    git clone https://github.com/mdolab/romfoam && \
    cd romfoam && \
    ./Allwmake

#. **Source code**

   The following assumes you want to compile the ROMFoam packages from scratch. If you use the pre-compiled version mentioned above, skip this.
   
   ROMFoam depends on multiple prerequisites, to install ROMFoam

   Install MACH framework following this website https://dafoam.readthedocs.io/en/latest/Installation.html **NOTE**: please follow the folder structure suggested in the installation documentation. Install all the dependencies except for DAFoam (no need to run its regression test either). This is because we have included DAFoam in ROMFoam. So DAFoam will be compiled along with ROMFoam later.
 
   Download Slepc-3.11.2 at http://slepc.upv.es/download/distrib/slepc-3.11.2.tar.gz  Move the downloaded tarball to your $HOME/packages and untar it. Then go into the slepc-3.11.2 folder and run::

    ./configure && make SLEPC_DIR=$PWD

   Append the following to your $HOME/.bashrc file::
  
    # SLEPC
    export SLEPC_DIR=$HOME/packages/slepc-3.11.2
    export SLEPC_ARCH=real-opt
    export PATH=$SLEPC_DIR/$SLEPC_ARCH/bin:$PATH
    export PATH=$SLEPC_DIR/$SLEPC_ARCH/include:$PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SLEPC_DIR/$SLEPC_ARCH/lib
    export SLEPC_LIB=$SLEPC_DIR/$SLEPC_ARCH/lib

   Then source it::

    source $HOME/.bashrc
 
   Now, download the ROMFoam repo, rename it to romfoam and put it in the $HOME/repos folder, load the OpenFOAM environment::
 
    source $HOME/OpenFOAM/OpenFOAM-v1812/etc/bashrc

   and go into the $HOME/repos/romfoam folder, and run::
   
    ./Allwmake

   Finally, verify the installation by going to romfoam/python/reg_test/AhmedBody/runROM, and run::

    . ./Allclean.sh && sh ./run_reg_tests.sh

   The regression test should take about 10 minutes and you should see::

    Regression test: Linear ROM...
    Success!
    Regression test: Nonlinear ROM...
    Success!

   If the test takes more than 10 minutes or any of the test fails, check the log file in reg_file_* for details. Make sure you pass the regression test before using ROMFoam!