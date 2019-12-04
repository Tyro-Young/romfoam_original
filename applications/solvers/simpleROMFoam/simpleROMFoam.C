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
    
        }
        
        Info<< "End\n" << endl;
    
        adjObj.writeObjFuncValues();
    
        adjDev->calcFlowResidualStatistics("print",1); // 1: write residual variables
    }
    else 
    {
        Info<<"-mode "<<mode<<" arg not valid!"<<endl
            <<"Options are:"<<endl
            <<"flow"<<endl
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
