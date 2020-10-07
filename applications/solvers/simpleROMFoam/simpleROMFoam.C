/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v1.0

    Description:
    ROM for simpleFoam

\*---------------------------------------------------------------------------*/

static char help[] = "Solves a linear system in parallel with KSP in OpenFOAM.\n\n";

#include <slepcsvd.h>
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "interpolationCellPoint.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "AdjointIO.H"
#include "AdjointSolverRegistry.H"
#include "AdjointRASModel.H"
#include "AdjointIndexing.H"
#include "AdjointJacobianConnectivity.H"
#include "AdjointObjectiveFunction.H"
#include "AdjointDerivative.H"
#include "ReducedOrderModeling.H"
#include "nearWallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    argList::addOption
    (
        "mode",
        "flow",
        "Which mode to run"
        // options: flow, offlineLinear, onlineLinear, offlineNonlinear, or onlineNonlinear
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"

    // Initialize the petsc solver. This needs to be called after the case
    // setup so that petsc uses the OpenFOAM MPI_COMM
    SlepcInitialize(&argc,&argv,(char*)0,help);

    // read options
    word mode="flow";
    if (args.optionFound("mode"))
    {
        mode = word(args.optionLookup("mode")());
    }
    else
    {
        Info<<"mode arg not found! Use the default -mode flow"<<endl
            <<"Options are:"<<endl
            <<"flow"<<endl
            <<"evalObj"<<endl
            <<"offlineLinear"<<endl
            <<"onlineLinear"<<endl
            <<"offlineNonlinear"<<endl
            <<"onlineNonlinear"<<endl
            <<"Example: simpleROMFoam -mode flow"<<endl;
    }

    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    if (mode=="offlineLinear")
    {
        ReducedOrderModeling rom(runTime,mesh,adjIO,adjReg(),adjRAS(),adjIdx,adjCon(),adjObj,adjDev(),"linear");
        rom.solveOfflineLinear();
    }
    else if (mode=="onlineLinear")
    {
        ReducedOrderModeling rom(runTime,mesh,adjIO,adjReg(),adjRAS(),adjIdx,adjCon(),adjObj,adjDev(),"linear");
        rom.solveOnlineLinear();
    }
    else if (mode=="offlineNonlinear")
    {
        ReducedOrderModeling rom(runTime,mesh,adjIO,adjReg(),adjRAS(),adjIdx,adjCon(),adjObj,adjDev(),"nonlinear");
        rom.solveOfflineNonlinear();
    }
    else if (mode=="onlineNonlinear")
    {
        ReducedOrderModeling rom(runTime,mesh,adjIO,adjReg(),adjRAS(),adjIdx,adjCon(),adjObj,adjDev(),"nonlinear");
        rom.solveOnlineNonlinear();
    }
    else if (mode=="flow")
    {
        turbulence->validate();

        Info<< "\nStarting time loop\n" << endl;

        while (simple.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            // --- Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }

            laminarTransport.correct();
            turbulence->correct();

            adjObj.printObjFuncValues();

            runTime.write();

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            if(runTime.writeTime())
            {
                // 1: write residual variables
                adjDev->copyStates("Var2Ref");
                adjDev->calcFlowResidualStatistics("print",1);
            }

        }

        Info<< "End\n" << endl;

        adjObj.writeObjFuncValues();


        // use desired patch
        const polyPatch& body = mesh.boundaryMesh()[6];
        scalarList x, y, z; //stores list of local coordinates


        forAll(body.localPoints(), i)
        {
          scalar yVal = body.localPoints()[i].y();
          if (yVal > 0.09 & yVal < 0.11) // find points around desired area, need to know coordinates
          {
            x.append(body.localPoints()[i].x());
            y.append(body.localPoints()[i].y());
            z.append(body.localPoints()[i].z());
          }
        }

        //if running in parallel
        if (Pstream::parRun())
        {
          //find minimum indices from points found
          scalar indX = findMin(x);
          scalar indY = findMin(y);
          scalar indZ = findMin(z);

          scalar minX, minY, minZ;


          //this function returns -1 if min not found. this will happen if the selected
          //patch does not have points on a processor
          if (indX != -1)
          {
            minX = x[indX];
          }

          if (indY != -1)
          {
            minY = y[indY];
          }

          if (indZ != -1)
          {
            minZ = z[indZ];
          }

          //find the global min of each coordinate
          reduce(minX,minOp<scalar>());
          reduce(minY,minOp<scalar>());
          reduce(minZ,minOp<scalar>());


          label globalSize = x.size();
          label localSize = x.size(); //size of coordinate list on each proc

          reduce(globalSize,sumOp<label>()); //find total number of coordinates

          scalarList globalX, globalY, globalZ; //stpres global coordinates

          globalX.resize(globalSize,minX); //create lists initialized to min vals
          globalY.resize(globalSize,minY);
          globalZ.resize(globalSize,minZ);

          labelList procList; //stores what processor each coordinate is on
          procList.resize(globalSize,0);

          label nProcs = Pstream::nProcs();

          labelList startIdx, nLocalPoints;
          startIdx.resize(nProcs,0); //stores starting index when constructing global coordinate list
          nLocalPoints.resize(nProcs,0); //stores how many points are on each processor

          forAll(nLocalPoints,i)
          {
            if (Pstream::myProcNo() == i)
            {
              nLocalPoints[i] = localSize;
            }
          }

          Pstream::listCombineGather(nLocalPoints, maxEqOp<label>());
          Pstream::listCombineScatter(nLocalPoints);

          for (label i = 1; i < nProcs; i++)
          {
            if (Pstream::myProcNo() == i)
            {
              startIdx[i] = startIdx[i-1] + nLocalPoints[i-1];
            }
            Pstream::listCombineGather(startIdx, maxEqOp<label>());
            Pstream::listCombineScatter(startIdx);
          }

          label procNo = Pstream::myProcNo();
          label idx = startIdx[procNo]; //find correct starting index for proc

          if (x.size() != 0) //if list is not empty locally
          {
            forAll(x,i)
            {
              globalX[idx+i] = x[i];
              globalY[idx+i] = y[i];
              globalZ[idx+i] = z[i];
              procList[idx+i] = procNo;
            }
          }

          Pstream::listCombineGather(procList, maxEqOp<label>());
          Pstream::listCombineScatter(procList);

          Pstream::listCombineGather(globalX, maxEqOp<scalar>());
          Pstream::listCombineScatter(globalX);

          Pstream::listCombineGather(globalY, maxEqOp<scalar>());
          Pstream::listCombineScatter(globalY);

          Pstream::listCombineGather(globalZ, maxEqOp<scalar>());
          Pstream::listCombineScatter(globalZ);

          //now all global lists have been created

          scalar maxX = globalX[findMax(globalX)];

          label nInt = 10; //number of intervals you want to find a probe point on

          scalar dx = (maxX-minX)/nInt;

          scalarList zVals; //height coordinates of candidate probe points in interval
          labelList iVals, idxList; //indices of candidates, indices of desired points

          idxList.resize(nInt);

          for (label i = 0; i < nInt; i++)
          {
            scalar lowX = minX + i*dx;
            scalar highX = minX + (i+1)*dx;

            for (label k = 0; k < globalX.size(); k++)
            {
              scalar xk = globalX[k];

              if (lowX <= xk && highX > xk)
              {
                zVals.append(globalZ[k]);
                iVals.append(k);
              }
            }
            idxList[i] = iVals[findMax(zVals)]; //keep highest z value (along top of body)
            //clear candidate lists before next iteration
            zVals.clear();
            iVals.clear();
          }

          List<vector> probePoints; //list of coords of probe points
          labelList probeProcs; //list of procs of coord points
          probePoints.resize(nInt);
          probeProcs.resize(nInt);

          forAll(idxList,i) //for all found probe points
          {
            label idx = idxList[i];
            vector pointI(globalX[idx],globalY[idx],globalZ[idx]);
            probePoints[i] = pointI;
            probeProcs[i] = procList[idx];
          }

          interpolationCellPoint<scalar> interp(p); //using cell point interpolation on pressure field

          std::ofstream myFile;
          myFile.open("pressure.txt");

          forAll(probePoints,i)
          {
            vector pointI = probePoints[i];
            label procNo = probeProcs[i];
            scalar pInterp = 0;
            if (Pstream::myProcNo() == procNo)
            {
              label cellI = mesh.findCell(pointI);
              pInterp = interp.interpolate(pointI,cellI); //probed pressure value
            }
            reduce(pInterp,sumOp<scalar>());
            myFile << pInterp << std::endl;
          }

        }
        //if running in serial
        else
        {
          scalar minX = x[findMin(x)];
          scalar maxX = x[findMax(x)];

          label nInt = 10;

          scalar dx = (maxX-minX)/nInt;

          scalarList zVals;
          labelList iVals, idxList;

          for (label i = 0; i < nInt; i++)
          {
            scalar lowX = minX + i*dx;
            scalar highX = minX + (i+1)*dx;

            for (label k = 0; k < x.size(); k++)
            {
              scalar xk = x[k];

              if (lowX <= xk && highX > xk)
              {
                zVals.append(z[k]);
                iVals.append(k);
              }
            }
            idxList.append(iVals[findMax(zVals)]);
            zVals.clear();
            iVals.clear();
          }

          List<vector> probePoints;
          probePoints.resize(nInt);

          forAll(idxList,i)
          {
            label idx = idxList[i];
            vector pointI(x[idx],0,z[idx]);
            probePoints[i] = pointI;
          }

          interpolationCellPoint<scalar> interp(p);
          std::ofstream myFile;
          myFile.open("pressure.txt");

          forAll(probePoints,i)
          {
            vector pointI = probePoints[i];
            label cellI = mesh.findCell(pointI);
            myFile << interp.interpolate(pointI,cellI);
            myFile << std::endl;
          }
        }


    }
    else if (mode=="evalObj")
    {
        Info<<"Evaluating the objective functions..."<<endl;
        adjObj.printObjFuncValues();
        Info<<"Writting the objective functions to disk..."<<endl;
        adjObj.writeObjFuncValues();
    }
    else
    {
        Info<<"-mode "<<mode<<" arg not valid!"<<endl
            <<"Options are:"<<endl
            <<"flow"<<endl
            <<"evalObj"<<endl
            <<"offlineLinear"<<endl
            <<"onlineLinear"<<endl
            <<"offlineNonlinear"<<endl
            <<"onlineNonlinear"<<endl
            <<"Example: simpleROMFoam -mode flow"<<endl;
    }


    SlepcFinalize();

    return 0;
}


// ************************************************************************* //
