.. _Tutorials:

Tutorials
==========

Linear ROM
----------

Case summary
~~~~~~~~~~~~

.. list-table:: Summary
   :widths: 25 25 
   :header-rows: 1

   * - ROM solver
     - simpleROMFoam
   * - Geometry
     - Ahmed body
   * - Mesh
     - 12,800 cells 
   * - Function
     - CD
   * - Re
     - 1.3e6
   * - alpha training
     - [5, 10, 15, 25, 20] degree
   * - alpha training
     - [15, 21] degree
   * - Turbulence Model 
     - Spalart-Allmaras
   * - CPU cores
     - 2

How to run
~~~~~~~~~~

To run this case, go to romfoam/tutorial/AhmedBody/runROM and run::

   ./AllrunLinear.sh
  
A few notes:

- The default setup has one design variable (ramp angle). 
- It will first run the OpenFOAM's built in incompressible flow solver (simpleFoam) to generate samples with various ramp angle values (stored at /tutorial/AhmedBody/sample*). 
- Then it will run simpleROMFoam to perform SVD and compute the reduced dR/dW (Ar) and dR/dX (Br) matrices and save them to the disk. 
- Finally, it will run simpleROMFoam to predict new fields given a new set of design variables (stored at /tutorial/AhmedBody/prediction*). 
- We use the last sample :math:`\alpha=20^{\circ}` as the reference point

After the run, the linear ROM has the following folder structures

.. code-block:: bash

   AhmedBody
   |--runROM
   |--output
   |--sample1 # folders containing the CFD samples
   |--sample2
   |--sample3
   |--sample4
   |--sample5
   |--prediction1 # folders for the predictions
   |--prediction2

How to post-process
~~~~~~~~~~~~~~~~~~~

To post-process this case:

- Go to tutorial/AhmedBody/prediction1. 
- Open ParaView and load the paraview.foam file. 
- The flow field variables computed by simpleFoam have names such as "U", "p", etc., while the flow fields computed from the reduced-order modeling have names such as "UROM", "pROM", etc. The modes from the SVD have names such as "UMode0", "UMode1", etc.

Comparison between the CFD (left) and ROM (right) results

 .. image:: images/linearROM_U.png
  :width: 400
  :alt: Alternative text

 .. image:: images/linearROM_UROM.png
  :width: 400
  :alt: Alternative text

Velocity modes: mode 0 (left) and mode 1 (right)

 .. image:: images/linearROM_UMode0.png
  :width: 400
  :alt: Alternative text

 .. image:: images/linearROM_UMode1.png
  :width: 400
  :alt: Alternative text

Details of the ROM log file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- When the ROM runs, it will print some information on screen. The first portion of the log file is for mesh generation. Then, it will run the CFD samples from 5 to 25 degrees, the following is the output for running CFD for 5 degree ramp angle. The computed CD is 0.4033, and the flow runtime is 23 s

.. code-block:: bash

   +--------------------------------------------------------------------------+
   |                    Evaluating Objective Function                         |
   +--------------------------------------------------------------------------+
   ('Design Variables: ', {'rampAngle': array([ 5.+0.j])}) # 5 degree ramp
   Deleting Previous Solution Files
   DVGeo PointSet UpToDate: False
   Updating DVGeo PointSet.... # deform the surface mesh
   DVGeo PointSet UpToDate: True
   Warping the volume mesh.... # deform the volume mesh
   Writting the updated volume mesh....
   Checking Mesh Quality. # check the deformed mesh quality
   Copying checkMeshLog to ../output/checkMeshLog_000
   Checking Mesh Quality. Passed!
   Calling Flow Solver 000 # run CFD, the log is recorded in the flowLog file
   Simulation Started. Check the flowLog file for the progress.
   Simulation Finished!
   CD std: 8.97089522285e-10
   ('Objective Functions: ', {'fail': False, 'CD': 0.403322747714277}) # computed CD
   ('Flow Runtime: ', 23.00730609893799) # Runtime

- Next, it will writ the dXv/dX matrix deltaVolPointMatPlusEps.bin, this matrix will be used when computing the dRdFFD matrix and eventually the Br matrix

.. code-block:: bash

   write deltaVolPointsMat at sample=4
   {'rampAngle': array([ 20.+0.j])}
   Delete Old DeltaVolPointMat ... 
   epsFFD: 0.001
   Writting deltaVolPointMatPlusEps for rampAngle
   (4, ['rampAngle'])

- After this, it will deform the geometry to the reference sample point (20 degree in this case)

.. code-block:: bash

   Deforming at sample=4
   {'rampAngle': array([ 20.+0.j])}

- Now, it runs the offline ROM

.. code-block:: bash

   Setting the w snapshot matrix. 1 s # assemble the snapshot matrix
   Reading variables from sample 1
   Reading variables from sample 2
   Reading variables from sample 3
   Reading variables from sample 4
   Reading variables from sample 5
   Solving the SVD...  # compute the SVD
   Singular Values    
   --------------------- 
   10670.62189947045
   1146.282545401202
   144.7740615824653
   38.11703859539994
   10.40819727165931
   Writing the phiMat # save the Phi matrix
   Calculating dRdW and dRdFFD at time = 5
   Computing dRdW*Phi...
   Computing phiT*dRdW*phi # compute Ar
   Calculating dRdFFD... 
   Initializing the dRdFFD matrix. 1 s
   dRdFFD matrix Created. 1 s
   Reading deltaVolPointMat
   Calculating dRdFFD...
   Computing phiT*dRdFFD # compute Br
   Writing the reduced matrices....  # save Ar and Br
   Writing the reduced matrices.... Done!

- Right after the offline ROM we need to clean up everything and run CFD for the reference sample, this is to provide the reference values for the onlineROM

.. code-block:: bash

   +--------------------------------------------------------------------------+
   |                    Evaluating Objective Function                         |
   +--------------------------------------------------------------------------+
   ('Design Variables: ', {'rampAngle': array([ 20.+0.j])})
   Deleting Previous Solution Files
   DVGeo PointSet UpToDate: False
   Updating DVGeo PointSet....
   DVGeo PointSet UpToDate: True
   Warping the volume mesh....
   Writting the updated volume mesh....
   Checking Mesh Quality.
   Copying checkMeshLog to ../output/checkMeshLog_000
   Checking Mesh Quality. Passed!
   Calling Flow Solver 000
   Simulation Started. Check the flowLog file for the progress.
   Simulation Finished!
   Copying flowLog to ../output/flowLog_000
   Copying objFuncs.dat to ../output/objFuncs_000.dat
   Reading ../output/objFuncs_000.dat
   ('keys', ['CD'])
   CD std: 5.93870308882e-08
   ('Objective Functions: ', {'fail': False, 'CD': 0.4009651757371485})
   ('Flow Runtime: ', 28.162258863449097)

- Finally, it runs the online ROM

.. code-block:: bash

   # first read the Phi, Ar, and Br matrices
   Initializing dRdWReduced matrix. 0 s
   dRdWReduced matrix Created. 0 s
   Initializing the dRdFFDReduced matrix. 0 s
   dRdFFDReduced matrix Created. 0 s
   Initializing the Phi matrix for W. 0 s
   Phi matrix Created. 0 s
   Solver Type: gmres
   GMRES Restart: 500
   ASM Overlap: 1
   Global PC Iters: 0
   Local PC Iters: 1
   Mat ReOrdering: rcm
   ILU PC Fill Level: 0
   GMRES Max Iterations: 500
   GMRES Relative Tolerance: 1e-06
   GMRES Absolute Tolerance: 1e-16 # solve the online reduced equation
   Main iteration 0 KSP Residual norm 3.662424193474e+02 0 s
   Main iteration 1 KSP Residual norm 5.312838779836e-14 0 s
   Main iteration 1 KSP Residual norm 5.312838779836e-14 0 s 
   Total iterations 1
   Writing Objective Function Values to objFuncs.dat
   CD: 0.3992995488257811 # CD from RoM
   Run the flow at sample = 1
   CD 0.3948548548033454  # CD from CFD


Details of the AllrunLinear.sh script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users need to prescribe the following parameters in the AllrunLinear.sh script, and they usually don't need to touch the rest

.. code-block:: bash

   #!/bin/bash
   exec=mpirun # mpi executive
   nProcs=2 # number of CPU cores
   solver=simpleROMFoam # solver name
   runEndTime=500 # how many steps to run for CFD
   nSamples=5 # number of samples
   refSample=$nSamples # use which reference sample, always set it to the last sample
   nDVs=1 # number of design variables
   predictSamples="1 2" # which indices to predict in the onlineROM, the details are defined by the DVs_Predict list in runFlow.py

Details of the runFlow.py script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are many parameters in the runFlow.py script. But the DVs_Train and DVs_Predict lists are the two users usually want to change. They define the values of the design variables for CFD samples and ROM prediction, respectively.

.. code-block:: bash

   if optVars[0]=='rampAngle':
   # 5 samples, 1 DV
   # we have five samples, their rampAngle is from 5 to 25 degrees, the sequence of the CFD samples does not matter but always put the reference sample as the last element in the DVs_Train list
   DVs_Train=   [[5.0],
                 [10.0],
                 [15.0],
                 [25.0],
                 [20.0]]
   # we want to predict two points at 15 and 21 degrees, this needs to be consistent with the predictSamples="1 2" parameter in AllrunLinear.sh script
   DVs_Predict= [[15.0], 
                 [21.0]]

Modifications: How to increase the number of CFD samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To increase the CFD samples from 5 to 6, we need to 

- Change nSamples=5 to nSamples=6 in AllrunLinear.sh
- In runFlow.py, append [22.0] to the DVs\_Train list

.. code-block:: bash

   if optVars[0]=='rampAngle':
     # 6 samples, 1 DV
     DVs_Train=   [[5.0],
                   [10.0],
                   [15.0],
                   [25.0],
                   [22.0],
                   [20.0]]
     DVs_Predict= [[15.0],
                   [21.0]]
  NOTE: we still put [20.0] as the last sample because we want it to be the reference sample

- Run ``./Allclean.sh`` and then run ``./AllrunLinear.sh``

Modifications: How to increase the number of predictions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To increase the predictions from 2 to 3, we need to 

- Change predictSamples="1 2" to predictSamples="1 2 3" in AllrunLinear.sh
- In runFlow.py, append [23.0] to the DVs\_Predict list

.. code-block:: bash

   if optVars[0]=='rampAngle':
     # 6 samples, 1 DV
     DVs_Train=   [[5.0],
                   [10.0],
                   [15.0],
                   [25.0],
                   [20.0]]
     DVs_Predict= [[15.0],
                   [21.0],
                   [23.0]]

- Run ``./Allclean.sh`` and then run ``./AllrunLinear.sh``

Modifications: How to use more CPU cores and run flow for more steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use 4 CPU cores and run the flow for 1000 steps, we need to 

- Change nProcs=2 to nProcs=4 in AllrunLinear.sh
- Change runEndTime=500 to runEndTime=1000 in AllrunLinear.sh
- Run ``./Allclean.sh`` and then run ``./AllrunLinear.sh``

Modifications: How to use more mesh cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To increase the number of mesh cells, we need to

- In system/blockMeshDict, increase the density of background mesh by changing::

   \texttt{hex (0 1 2 3 4 5 6 7) (20 5 5) simpleGrading (1 1 1)} 

  to::     

   \texttt{hex (0 1 2 3 4 5 6 7) (40 10 10) simpleGrading (1 1 1)}

- In system/snappyHexMeshDict, add prism boundary layer mesh by setting ``addLayers true;``, and increase the desired number of prism layers by setting ``nSurfaceLayers`` to a higher number

- In system/snappyHexMeshDict, you can also increase the mesh refinement levels in castellatedMeshControls-features and castellatedMeshControls-refinementSurfaces from 4 to 5

- Check this website for details on how to generate mesh using snappyHexMesh: https://cfd.direct/openfoam/user-guide/v6-snappyhexmesh


Modifications: How to use more design variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run ROM with two design variables (ride height and ramp angle case), we need to

- Change nDVs=1 to nDVs=2, and nSamples=5 to nSamples=10 in AllrunLinear.sh. Note: in runFlow.py we have 10 samples, see following

.. code-block:: bash

   elif optVars[0]=='rampAngleAndRideHeight':
   # 10 samples, 2 DVs
   DVs_Train=   [[24.1,0.045],
                 [25.3,0.049],
                 [25.7,0.043],
                 [25.1,0.059],
                 [25.5,0.055],
                 [24.5,0.041],
                 [24.3,0.053],
                 [24.7,0.057],
                 [25.9,0.047],
                 [24.9,0.051]]
   DVs_Predict= [[25.7,0.043],
                 [24.1,0.045]]

- In runFlow.py, change the default value from ``default="['rampAngle']"`` to ``default="['rampAngleAndRideHeight']"`` for the ``--optVars`` option 
- Run ``./Allclean.sh`` and then run ``./AllrunLinear.sh``

Nonlinear ROM
-------------

How to run
~~~~~~~~~~~

We use the same tutorial as the linear ROM case. 
  
To use nonlinear ROM, go to romfoam/tutorial/AhmedBody/runROM and::

   ./AllrunNonlinear.sh

  
A few notes:
 
- The overall process is similar to the linear ROM
- In the offline stage, we no longer need to compute Ar and Br, instead we need to compute the preconditioner matrix :math:`\frac{\partial \overrightarrow{R}_r}{\partial \overrightarrow{w}_r}_{\textrm{PC}}` and save it to the disk
- The nonlinear ROM is implemented in a brute-force manner and the online stage is generally slower than that in the linear ROM. In addition, we found that the convergence of the Newton method depends on the initial flow fields. Generally, we should use a sample point that is closest to the predicted point as initial fields

Details of the ROM log file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mesh generation and offline ROM logs are same as the linear ROM. For the onlineROM, the Newton-Krylov solution output is

.. code-block:: bash

   -----------------------------------------------------------------------------------
   |Main |Func | Step |  Lin  | RunTime |  Res Norm|  Res Norm|  Res Norm |     CD 
   |Iter |Eval |      |  Res  | (s)     |  (Turb)  |  (Full)  |  (Reduced)|           
   -----------------------------------------------------------------------------------
       0     0 1.0000 0.0e+00 1.0000e+00 1.0675e+00 1.8301e+03 3.6981e+02  3.882030e-01
       1     5 0.5000 6.1e-04 1.0000e+00 1.5580e+00 1.0280e+03 2.2030e+02  3.924965e-01
       2    10 1.0000 2.4e-04 1.0000e+00 1.6684e+00 6.9184e+02 2.0213e+02  3.962898e-01
       3    14 1.0000 5.2e-03 1.0000e+00 4.1845e-02 1.4301e+01 4.4271e+00  3.948602e-01
       4    19 1.0000 8.0e-04 1.0000e+00 3.2940e-04 2.4071e-01 5.9697e-03  3.948558e-01
       5    25 1.0000 0.0e+00 2.0000e+00 3.0532e-04 1.7259e-01 5.5937e-08  3.948548e-01
       6    31 1.0000 0.0e+00 2.0000e+00 3.0532e-04 1.7259e-01 4.0091e-12  3.948548e-01
   Absolute Tolerance 4.009128232387992e-12 less than the presribed nkAbsTol 1e-08
   NK completed!
   Writing Objective Function Values to objFuncs.dat
   CD: 0.3948548668384586 # CD predict by ROM
   Run the flow at sample = 1
   CD 0.394854972380792  # CD predict by CFD

**NOTE:** Res Norm (Reduced) is the reduced residual and Res Norm (Full) is the full-scale residual, here we drive only the reduced residual, and with the hope that the full-scale residual also drops (it does drop in this case).

  