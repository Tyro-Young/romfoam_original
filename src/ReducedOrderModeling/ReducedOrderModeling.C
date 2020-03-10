/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    ROM stuff

\*---------------------------------------------------------------------------*/

#include "ReducedOrderModeling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

ReducedOrderModeling::ReducedOrderModeling
(
    Foam::Time& runTime,
    fvMesh& mesh,
    AdjointIO& adjIO,
    AdjointSolverRegistry& adjReg,
    AdjointRASModel& adjRAS,
    AdjointIndexing& adjIdx,
    AdjointJacobianConnectivity& adjCon,
    AdjointObjectiveFunction& adjObj,
    AdjointDerivative& adjDev,
    word mode
)
    :
    runTime_(runTime),
    mesh_(mesh),
    adjIO_(adjIO),
    adjReg_(adjReg),
    adjRAS_(adjRAS),
    adjIdx_(adjIdx),
    adjCon_(adjCon),
    adjObj_(adjObj),
    adjDev_(adjDev),
    mode_(mode),
    db_(mesh.thisDb()),
    // read dict from system/romDict
    romDict_
    (
        IOobject
        (
            "romDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    romParameters_
    (
        IOobject
        (
            "romParameters",
            mesh_.time().system(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    )
{

    // read
    nSamples                  = readOptionOrDefault<label>(romDict_,"nSamples",1);
    deltaFFD                  = readOptionOrDefault<scalarList>(romDict_,"deltaFFD",{});
    svdType                   = readOptionOrDefault<word>(romDict_,"svdType","cross");
    svdTol                    = readOptionOrDefault<scalar>(romDict_,"svdTol",1e-8);
    svdMaxIts                 = readOptionOrDefault<label>(romDict_,"svdMaxIts",100);
    svdRequestedN             = readOptionOrDefault<label>(romDict_,"svdRequestedN",1);
    useMF                     = readOptionOrDefault<label>(romDict_,"useMF",1);
    debugMode                 = readOptionOrDefault<label>(romDict_,"debugMode",0);
    mfStep                    = readOptionOrDefault<scalar>(romDict_,"mfStep",1e-6);
    romNKAbsTol               = readOptionOrDefault<scalar>(romDict_,"romNKAbsTol",1e-8);
    romNKGMRESRTol            = readOptionOrDefault<scalar>(romDict_,"romNKGMRESRTol",1e-2);
    romNKGMRESMaxLS           = readOptionOrDefault<label>(romDict_,"romNKGMRESMaxLS",10);
    romNKMaxIts               = readOptionOrDefault<label>(romDict_,"romNKMaxIts",20);
    useLSPG                   = readOptionOrDefault<label>(romDict_,"useLSPG",0);
    romNKMFFDH                = readOptionOrDefault<scalar>(romDict_,"romNKMFFDH",-9999.0);
    

    // print all the parameters to screen    
    Info<<"ROM Parameters"<<romParameters_<<endl;

    localSize_=adjIdx_.nLocalAdjointStates;
    nFFDs_=adjIO_.nFFDPoints;

    label checkMode=0;
    if (mode_=="nonlinear") checkMode=1;
    if (mode_=="linear") checkMode=1;
    if(checkMode!=1) FatalErrorIn("")<<"mode should be linear or nonlinear"<< abort(FatalError);

}

ReducedOrderModeling::~ReducedOrderModeling()
{
}

void ReducedOrderModeling::initializeOnlineLinear()
{
    
    this->initializedRdWMatReduced();
    this->initializedRdFFDMatReduced();
    this->initializeSVDPhiMat();

    label nProcs = Pstream::nProcs();
    std::ostringstream np("");
    np<<nProcs;
    std::string fNamedRdWReduced="dRdWReduced_"+np.str();
    std::string fNamedRdFFDReduced="dRdFFDReduced_"+np.str();
    std::string fNamePhi="svdPhiWMat_"+np.str();

    adjIO_.readMatrixBinary(dRdWReduced_,fNamedRdWReduced);
    adjIO_.readMatrixBinary(dRdFFDReduced_,fNamedRdFFDReduced);
    adjIO_.readMatrixBinary(svdPhiWMat_,fNamePhi);
    
    VecCreate(PETSC_COMM_WORLD,&deltaFFDVec_);
    VecSetSizes(deltaFFDVec_,PETSC_DECIDE,nFFDs_);
    VecSetFromOptions(deltaFFDVec_);
    VecZeroEntries(deltaFFDVec_);

    label Istart,Iend;
    VecGetOwnershipRange(deltaFFDVec_,&Istart,&Iend);
    
    for(label i=Istart;i<Iend;i++)
    {
        scalar val=deltaFFD[i];
        VecSetValue(deltaFFDVec_,i,val,INSERT_VALUES);
    }
    VecAssemblyBegin(deltaFFDVec_);
    VecAssemblyEnd(deltaFFDVec_);

    if(debugMode) adjIO_.writeVectorASCII(deltaFFDVec_,"deltaFFDVec");

    MatCreate(PETSC_COMM_WORLD,&dRdWPhi_);
    MatSetSizes(dRdWPhi_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(dRdWPhi_);
    MatMPIAIJSetPreallocation(dRdWPhi_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(dRdWPhi_,nSamples,NULL);
    MatSetUp(dRdWPhi_);

}

void ReducedOrderModeling::initializeOnlineNonlinear()
{
    
    this->initializeSVDPhiMat();

    // read svdPhiWMat
    label nProcs = Pstream::nProcs();
    std::ostringstream np("");
    np<<nProcs;
    std::string fNamePhi="svdPhiWMat_"+np.str();
    adjIO_.readMatrixBinary(svdPhiWMat_,fNamePhi);

    if(useLSPG==0)
    {
        // read svdPhiRMat
        std::ostringstream npR("");
        npR<<nProcs;
        std::string fNamePhiR="svdPhiRMat_"+npR.str();
        adjIO_.readMatrixBinary(svdPhiRMat_,fNamePhiR);
    }
    
    VecCreate(PETSC_COMM_WORLD,&deltaFFDVec_);
    VecSetSizes(deltaFFDVec_,PETSC_DECIDE,nFFDs_);
    VecSetFromOptions(deltaFFDVec_);
    VecZeroEntries(deltaFFDVec_);

    label Istart,Iend;
    VecGetOwnershipRange(deltaFFDVec_,&Istart,&Iend);
    
    for(label i=Istart;i<Iend;i++)
    {
        scalar val=deltaFFD[i];
        VecSetValue(deltaFFDVec_,i,val,INSERT_VALUES);
    }
    VecAssemblyBegin(deltaFFDVec_);
    VecAssemblyEnd(deltaFFDVec_);

    if(debugMode) adjIO_.writeVectorASCII(deltaFFDVec_,"deltaFFDVec");

    VecCreate(PETSC_COMM_WORLD,&wVecFull_);
    VecSetSizes(wVecFull_,localSize_,PETSC_DECIDE);
    VecSetFromOptions(wVecFull_);
    VecZeroEntries(wVecFull_);

    VecDuplicate(wVecFull_,&rVecFull_);
    VecZeroEntries(rVecFull_);

    VecCreate(PETSC_COMM_WORLD,&wVecReduced_);
    VecSetSizes(wVecReduced_,PETSC_DECIDE,nSamples);
    VecSetFromOptions(wVecReduced_);
    VecZeroEntries(wVecReduced_);

    VecDuplicate(wVecReduced_,&rVecReduced_);
    VecZeroEntries(rVecReduced_);

    MatCreate(PETSC_COMM_WORLD,&dRdWPhi_);
    MatSetSizes(dRdWPhi_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(dRdWPhi_);
    MatMPIAIJSetPreallocation(dRdWPhi_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(dRdWPhi_,nSamples,NULL);
    MatSetUp(dRdWPhi_);

}

void ReducedOrderModeling::initializeOfflineLinear()
{
    MatCreate(PETSC_COMM_WORLD,&dRdWPhi_);
    MatSetSizes(dRdWPhi_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(dRdWPhi_);
    MatMPIAIJSetPreallocation(dRdWPhi_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(dRdWPhi_,nSamples,NULL);
    MatSetUp(dRdWPhi_);
}

void ReducedOrderModeling::initializeOfflineNonlinear()
{

    VecCreate(PETSC_COMM_WORLD,&wVecFull_);
    VecSetSizes(wVecFull_,localSize_,PETSC_DECIDE);
    VecSetFromOptions(wVecFull_);
    VecZeroEntries(wVecFull_);

    VecDuplicate(wVecFull_,&rVecFull_);
    VecZeroEntries(rVecFull_);

    VecCreate(PETSC_COMM_WORLD,&wVecReduced_);
    VecSetSizes(wVecReduced_,PETSC_DECIDE,nSamples);
    VecSetFromOptions(wVecReduced_);
    VecZeroEntries(wVecReduced_);

    VecDuplicate(wVecReduced_,&rVecReduced_);
    VecZeroEntries(rVecReduced_);

    MatCreate(PETSC_COMM_WORLD,&dRdWPhi_);
    MatSetSizes(dRdWPhi_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(dRdWPhi_);
    MatMPIAIJSetPreallocation(dRdWPhi_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(dRdWPhi_,nSamples,NULL);
    MatSetUp(dRdWPhi_);

}

void ReducedOrderModeling::initializeSnapshotMat()
{
    Info<<"Initializing the Snapshot matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // now initialize the memory for the jacobian itself
    // create wSnapshotMat_
    MatCreate(PETSC_COMM_WORLD,&wSnapshotMat_);
    MatSetSizes(wSnapshotMat_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(wSnapshotMat_);
    MatMPIAIJSetPreallocation(wSnapshotMat_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(wSnapshotMat_,nSamples,NULL);
    MatSetUp(wSnapshotMat_);
    Info<<"W Snapshot matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;


    // create rSnapshotMat_
    if(mode_=="nonlinear" && useLSPG==0)
    {
        MatCreate(PETSC_COMM_WORLD,&rSnapshotMat_);
        MatSetSizes(rSnapshotMat_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
        MatSetFromOptions(rSnapshotMat_);
        MatMPIAIJSetPreallocation(rSnapshotMat_,nSamples,NULL,nSamples,NULL);
        MatSeqAIJSetPreallocation(rSnapshotMat_,nSamples,NULL);
        MatSetUp(rSnapshotMat_);
        Info<<"R Snapshot matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    }


    return;
}

void ReducedOrderModeling::initializedRdFFDMat()
{
    //******************************** Compute dRdFFD *****************************//
    Info<<"Initializing the dRdFFD matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // create dRdFFD_
    MatCreate(PETSC_COMM_WORLD,&dRdFFD_);
    MatSetSizes(dRdFFD_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nFFDs_);
    MatSetFromOptions(dRdFFD_);
    MatMPIAIJSetPreallocation(dRdFFD_,nFFDs_,NULL,nFFDs_,NULL);
    MatSeqAIJSetPreallocation(dRdFFD_,nFFDs_,NULL);
    MatSetOption(dRdFFD_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(dRdFFD_);
    Info<<"dRdFFD matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    return;
}

void ReducedOrderModeling::initializeSVDPhiMat()
{
    Info<<"Initializing the Phi matrix for W. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // create svdPhiWMat_
    MatCreate(PETSC_COMM_WORLD,&svdPhiWMat_);
    MatSetSizes(svdPhiWMat_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(svdPhiWMat_);
    MatMPIAIJSetPreallocation(svdPhiWMat_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(svdPhiWMat_,nSamples,NULL);
    MatSetOption(svdPhiWMat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(svdPhiWMat_);
    Info<<"Phi matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    if(mode_=="nonlinear" && useLSPG==0)
    {
        Info<<"Initializing the Phi matrix for R. "<<runTime_.elapsedClockTime()<<" s"<<endl;
        // create svdPhiRMat_
        MatCreate(PETSC_COMM_WORLD,&svdPhiRMat_);
        MatSetSizes(svdPhiRMat_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
        MatSetFromOptions(svdPhiRMat_);
        MatMPIAIJSetPreallocation(svdPhiRMat_,nSamples,NULL,nSamples,NULL);
        MatSeqAIJSetPreallocation(svdPhiRMat_,nSamples,NULL);
        MatSetOption(svdPhiRMat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        MatSetUp(svdPhiRMat_);
        Info<<"Phi matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    }

    return;
}

void ReducedOrderModeling::initializedRdWMatReduced()
{
    Info<<"Initializing dRdWReduced matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    MatCreate(PETSC_COMM_WORLD,&dRdWReduced_);
    MatSetSizes(dRdWReduced_,PETSC_DECIDE,PETSC_DECIDE,nSamples,nSamples);
    MatSetFromOptions(dRdWReduced_);
    MatMPIAIJSetPreallocation(dRdWReduced_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(dRdWReduced_,nSamples,NULL);
    MatSetOption(dRdWReduced_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(dRdWReduced_);

    Info<<"dRdWReduced matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    return;
}

void ReducedOrderModeling::initializedRdFFDMatReduced()
{
    Info<<"Initializing the dRdFFDReduced matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    MatCreate(PETSC_COMM_WORLD,&dRdFFDReduced_);
    MatSetSizes(dRdFFDReduced_,PETSC_DECIDE,PETSC_DECIDE,nSamples,nFFDs_);
    MatSetFromOptions(dRdFFDReduced_);
    MatMPIAIJSetPreallocation(dRdFFDReduced_,nFFDs_,NULL,nFFDs_,NULL);
    MatSeqAIJSetPreallocation(dRdFFDReduced_,nFFDs_,NULL);
    MatSetOption(dRdFFDReduced_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(dRdFFDReduced_);

    Info<<"dRdFFDReduced matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    return;
}

void ReducedOrderModeling::setSnapshotMat()
{
    // NOTE: make sure you copy the snap shot flow fields and rename them to 1, 2, 3, etc.
    // and put them in the main folder where you run offlineROM

    Info<<"Setting the w snapshot matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    MatZeroEntries(wSnapshotMat_);

    for(label idxI=0; idxI<nSamples;idxI++)
    {
        word varDir = name(idxI+1);
        Info<< "Reading variables from sample " <<varDir<< endl;

        forAll(adjReg_.volVectorStates,idxJ)                                           
        {        
            const word stateName = adjReg_.volVectorStates[idxJ];                      
            volVectorField state
            (
                IOobject
                (
                    stateName,
                    varDir,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false // not registering the IOobject
                ),
                mesh_
            );

            forAll(state,cellI)
            {
                for(label i=0;i<3;i++)
                {
                    PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex(stateName,cellI,i);
                    PetscInt colI = idxI;
                    PetscScalar val = state[cellI][i];
                    MatSetValues(wSnapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
                }
            }                                                           
        }
    
        forAll(adjReg_.volScalarStates,idxJ)
        {
            const word stateName = adjReg_.volScalarStates[idxJ];                      
            volScalarField state
            (
                IOobject
                (
                    stateName,
                    varDir,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false // not registering the IOobject
                ),
                mesh_
            );

            forAll(state,cellI)
            {
                PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex(stateName,cellI);
                PetscInt colI = idxI;
                PetscScalar val = state[cellI];
                MatSetValues(wSnapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
            }          
        }
    
        forAll(adjRAS_.turbStates,idxJ)
        {
            const word stateName = adjRAS_.turbStates[idxJ];                      
            volScalarField state
            (
                IOobject
                (
                    stateName,
                    varDir,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false // not registering the IOobject
                ),
                mesh_
            );

            forAll(state,cellI)
            {
                PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex(stateName,cellI);
                PetscInt colI = idxI;
                PetscScalar val = state[cellI];
                MatSetValues(wSnapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
            }          
    
        }

        forAll(adjReg_.surfaceScalarStates,idxJ)
        {
            const word stateName = adjReg_.surfaceScalarStates[idxJ];                      
            surfaceScalarField state
            (
                IOobject
                (
                    stateName,
                    varDir,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false // not registering the IOobject
                ),
                mesh_
            );

            forAll(mesh_.faces(), faceI)
            {
                PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex(stateName,faceI);
                PetscInt colI = idxI;
                PetscScalar val;
                if (faceI < adjIdx_.nLocalInternalFaces)
                {
                    val = state[faceI];
                }
                else
                {
                    label relIdx=faceI-adjIdx_.nLocalInternalFaces;
                    label patchIdx=adjIdx_.bFacePatchI[relIdx];
                    label faceIdx=adjIdx_.bFaceFaceI[relIdx];
                    val = state.boundaryField()[patchIdx][faceIdx];
                } 
                MatSetValues(wSnapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
            }
        }

    }

    MatAssemblyBegin(wSnapshotMat_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(wSnapshotMat_,MAT_FINAL_ASSEMBLY);

    if(debugMode)
    {
        adjIO_.writeMatrixASCII(wSnapshotMat_,"wSnapshotMat");
        adjIO_.writeMatrixBinary(wSnapshotMat_,"wSnapshotMat");
    }

    Info<<"The w snapshot matrix is set. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    if(mode_=="nonlinear"  && useLSPG==0)
    {
        Info<<"Setting the r snapshot matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;

        MatZeroEntries(rSnapshotMat_);
    
        for(label idxI=0; idxI<nSamples;idxI++)
        {
            word varDir = name(idxI+1);
            Info<< "Reading variables from sample " <<varDir<< endl;
    
            forAll(adjReg_.volVectorStates,idxJ)                                           
            {        
                const word stateName = adjReg_.volVectorStates[idxJ];
                const word stateResName = stateName+"Res";                      
                volVectorField stateRes
                (
                    IOobject
                    (
                        stateResName,
                        varDir,
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false // not registering the IOobject
                    ),
                    mesh_
                );
    
                forAll(stateRes,cellI)
                {
                    for(label i=0;i<3;i++)
                    {
                        PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex(stateName,cellI,i);
                        PetscInt colI = idxI;
                        PetscScalar val = stateRes[cellI][i];
                        MatSetValues(rSnapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
                    }
                }                                                           
            }
        
            forAll(adjReg_.volScalarStates,idxJ)
            {
                const word stateName = adjReg_.volScalarStates[idxJ];   
                const word stateResName = stateName+"Res";                   
                volScalarField stateRes
                (
                    IOobject
                    (
                        stateResName,
                        varDir,
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false // not registering the IOobject
                    ),
                    mesh_
                );
    
                forAll(stateRes,cellI)
                {
                    PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex(stateName,cellI);
                    PetscInt colI = idxI;
                    PetscScalar val = stateRes[cellI];
                    MatSetValues(rSnapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
                }          
            }
        
            forAll(adjRAS_.turbStates,idxJ)
            {
                const word stateName = adjRAS_.turbStates[idxJ];  
                const word stateResName = stateName+"Res";                    
                volScalarField stateRes
                (
                    IOobject
                    (
                        stateResName,
                        varDir,
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false // not registering the IOobject
                    ),
                    mesh_
                );
    
                forAll(stateRes,cellI)
                {
                    PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex(stateName,cellI);
                    PetscInt colI = idxI;
                    PetscScalar val = stateRes[cellI];
                    MatSetValues(rSnapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
                }          
        
            }
    
            forAll(adjReg_.surfaceScalarStates,idxJ)
            {
                const word stateName = adjReg_.surfaceScalarStates[idxJ]; 
                const word stateResName = stateName+"Res";                     
                surfaceScalarField stateRes
                (
                    IOobject
                    (
                        stateResName,
                        varDir,
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false // not registering the IOobject
                    ),
                    mesh_
                );
    
                forAll(mesh_.faces(), faceI)
                {
                    PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex(stateName,faceI);
                    PetscInt colI = idxI;
                    PetscScalar val;
                    if (faceI < adjIdx_.nLocalInternalFaces)
                    {
                        val = stateRes[faceI];
                    }
                    else
                    {
                        label relIdx=faceI-adjIdx_.nLocalInternalFaces;
                        label patchIdx=adjIdx_.bFacePatchI[relIdx];
                        label faceIdx=adjIdx_.bFaceFaceI[relIdx];
                        val = stateRes.boundaryField()[patchIdx][faceIdx];
                    } 
                    MatSetValues(rSnapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
                }
            }
    
        }
    
        MatAssemblyBegin(rSnapshotMat_,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(rSnapshotMat_,MAT_FINAL_ASSEMBLY);
    
        if(debugMode)
        {
            adjIO_.writeMatrixASCII(rSnapshotMat_,"rSnapshotMat");
            adjIO_.writeMatrixBinary(rSnapshotMat_,"rSnapshotMat");
        }
    
        Info<<"The r snapshot matrix is set. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    }
}

void ReducedOrderModeling::solveSVD(Mat snapshotMat, Mat phiMat, word phiMatName)
{
    PetscInt       nsv,maxit,nconv,its;
    PetscScalar    tol,sigma,val;
    SVDType        type;

    /// left vector for SVD uVec*PhiMat*vVec with size same as the row size for PhiMat
    Vec uVec;
    /// right vector for SVD, uVec*PhiMat*vVec with size same as the column size for PhiMat
    Vec vVec;
    // vVec has a size same as the colum size of snapshotMat, 
    // while uVec has a size same as the row size of snapshotMat
    MatCreateVecs(snapshotMat,&vVec,&uVec);
    VecZeroEntries(vVec);
    VecZeroEntries(uVec);

    // ----------- compute phiMat for W ----------------
    PetscInt Istart,Iend;
    MatGetOwnershipRange(phiMat,&Istart,&Iend);

    Info<<"Solving the SVD..."<<endl;
    SVD svd;
    SVDCreate(PETSC_COMM_WORLD,&svd);
    SVDSetOperator(svd,snapshotMat);
    SVDSetFromOptions(svd);

    if (svdType=="cross") SVDSetType(svd,SVDCROSS);
    else if (svdType=="cyclic") SVDSetType(svd,SVDCYCLIC);
    else if (svdType=="lapack") SVDSetType(svd,SVDLAPACK);
    else if (svdType=="lanczos") SVDSetType(svd,SVDTRLANCZOS);
    else if (svdType=="trlanczos") SVDSetType(svd,SVDTRLANCZOS);
    else
    {
        FatalErrorIn("")<<"svdType not supported!"<<
         "options are: cross, cyclic, lapack, lanczos, trlanczos"<< abort(FatalError);
    }

    SVDSetTolerances(svd,svdTol,svdMaxIts);
    SVDSetDimensions(svd,svdRequestedN,2*svdRequestedN,PETSC_DEFAULT);
    SVDSetImplicitTranspose(svd,PETSC_TRUE);

    SVDSolve(svd);

    // Output SVD solution diagnostics
    SVDGetIterationNumber(svd,&its);
    PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);
    SVDGetType(svd,&type);
    PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n",type);
    SVDGetDimensions(svd,&nsv,NULL,NULL);
    PetscPrintf(PETSC_COMM_WORLD," Number of requested singular values: %D\n",nsv);
    SVDGetTolerances(svd,&tol,&maxit);
    PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);
    SVDGetConverged(svd,&nconv);
    PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate singular triplets: %D\n\n",nconv);


    // Display singular values and relative errors
    PetscPrintf(PETSC_COMM_WORLD,
     "     Singular Values    \n"
     "  --------------------- \n");
    for (label i=0;i<nconv;i++) 
    {
        // Get converged singular triplets: i-th singular value is stored in sigma
        SVDGetSingularTriplet(svd,i,&sigma,uVec,vVec);
        Info<<"  "<<sigma<<endl;
        PetscScalar *uVecArray;
        VecGetArray(uVec,&uVecArray);
        for(label j=Istart;j<Iend;j++)
        {
            label relIdx = j-Istart;
            val = uVecArray[relIdx];
            MatSetValues(phiMat,1,&j,1,&i,&val,INSERT_VALUES);
        }
        VecRestoreArray(uVec,&uVecArray);
    }
    MatAssemblyBegin(phiMat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(phiMat,MAT_FINAL_ASSEMBLY);

    Info<<"Writing the phiMat"<<endl;
    label nProcs = Pstream::nProcs();
    std::ostringstream np("");
    np<<nProcs;
    std::string fNamePhi=phiMatName+"_"+np.str();

    adjIO_.writeMatrixBinary(phiMat,fNamePhi);
    if(debugMode) adjIO_.writeMatrixASCII(phiMat,fNamePhi);

    SVDDestroy(&svd);
}

void ReducedOrderModeling::solveOfflineNonlinear()
{
    this->initializeOfflineNonlinear();

    // save the unperturb residual statistics as reference
    // we will verify against this reference to ensure consistent residuals after perturbing
    // and resetting states
    adjDev_.calcFlowResidualStatistics("set");

    this->initializeSnapshotMat();

    adjDev_.calcFlowResidualStatistics("print");

    this->setSnapshotMat();

    this->initializeSVDPhiMat();

    this->solveSVD(wSnapshotMat_,svdPhiWMat_,"svdPhiWMat");

    this->getPhiMatStateInfo(svdPhiWMat_);

    if(useLSPG==0) this->solveSVD(rSnapshotMat_,svdPhiRMat_,"svdPhiRMat");

    if (debugMode) this->getMatrixColNorms(svdPhiWMat_,"svdPhiWMat");

    // initialize and calculate rdRdWPC
    Mat rdRdWPC;
    this->initializeReducedJacobian(&rdRdWPC);
    // save the unperturb residual statistics as reference
    this->calcReducedJacobian(rdRdWPC); 
    if (debugMode) adjIO_.writeMatrixASCII(rdRdWPC,"rdRdWPC");
    adjIO_.writeMatrixBinary(rdRdWPC,"rdRdWPC");

    Info<<"Done! "<<runTime_.elapsedClockTime()<<" s"<<endl;
}

void ReducedOrderModeling::solveOfflineLinear()
{

    this->initializeOfflineLinear();

    // save the unperturb residual statistics as reference
    // we will verify against this reference to ensure consistent residuals after perturbing
    // and resetting states
    adjDev_.calcFlowResidualStatistics("set");

    this->initializeSnapshotMat();

    adjDev_.calcFlowResidualStatistics("print");

    this->setSnapshotMat();

    this->initializeSVDPhiMat();

    this->solveSVD(wSnapshotMat_,svdPhiWMat_,"svdPhiWMat");

    if (debugMode) this->getMatrixColNorms(svdPhiWMat_,"svdPhiWMat");

    Info<<"Calculating dRdW and dRdFFD at time = "<<runTime_.elapsedClockTime()<<" s"<<endl;

    label nProcs = Pstream::nProcs();
    std::ostringstream np("");
    np<<nProcs;
    
    if(useMF) // matrix free
    {
        this->calcReducedMatsMF();
    }
    else
    {
        //******************************** Compute dRdW *****************************//
        // compute preallocation vecs for dRdW 
        adjCon_.setupdRdWCon(1);
        adjCon_.initializedRdWCon();
        adjCon_.setupdRdWCon(0);
        // Read in the precomputed coloring
        adjCon_.readdRdWColoring();
        
        label transposed=0,isPC=0;
        adjDev_.initializedRdW(&dRdW_,transposed);
    
        std::string fNamedRdW="dRdW_"+np.str();
        std::string fNamedRdWBin="dRdW_"+np.str()+".bin";
        std::ifstream fIn(fNamedRdWBin);
        if(fIn.fail()) 
        {
            Info<<"Calculating dRdW... "<<endl;
            adjDev_.calcdRdW(dRdW_,transposed,isPC);
            adjIO_.writeMatrixBinary(dRdW_,fNamedRdW);
        }
        else
        {
            Info<<"Reading dRdW... "<<endl;
            // read 
            PetscViewer viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD,fNamedRdWBin.c_str(),FILE_MODE_READ,&viewer);
            MatLoad(dRdW_,viewer);
            MatAssemblyBegin(dRdW_,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(dRdW_,MAT_FINAL_ASSEMBLY);
            PetscViewerDestroy(&viewer);
        }
    
        adjCon_.deletedRdWCon();
    
        //******************************** Compute dRdWReduced *****************************// 
        Info<< "Computing phiT*dRdW*phi" << endl;
        // now compute the matmat mult
        MatMatMult(dRdW_,svdPhiWMat_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdWPhi_);
        MatTransposeMatMult(svdPhiWMat_,dRdWPhi_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdWReduced_);
    
        MatDestroy(&dRdW_);
    
        this->initializedRdFFDMat();
        
        std::string fNamedRdFFD="dRdFFD_"+np.str();
        std::string fNamedRdFFDBin="dRdFFD_"+np.str()+".bin";
        std::ifstream fIn1(fNamedRdFFDBin);
        if(fIn1.fail()) 
        {
            Info<<"Calculating dRdFFD... "<<endl;
            adjDev_.calcdRdFFD(dRdFFD_);
            adjIO_.writeMatrixBinary(dRdFFD_,fNamedRdFFD);
        }
        else
        {
            Info<<"Reading dRdFFD... "<<endl;
            // read 
            PetscViewer viewer;
            PetscViewerBinaryOpen(PETSC_COMM_WORLD,fNamedRdFFDBin.c_str(),FILE_MODE_READ,&viewer);
            MatLoad(dRdFFD_,viewer);
            MatAssemblyBegin(dRdFFD_,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(dRdFFD_,MAT_FINAL_ASSEMBLY);
            PetscViewerDestroy(&viewer);
        }
        
        //******************************** Compute Br *****************************//
        Info<< "Computing phiT*dRdFFD" << endl;
        this->initializedRdFFDMatReduced();
        MatTransposeMatMult(svdPhiWMat_,dRdFFD_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdFFDReduced_);
    }

    Info<<"Writing the reduced matrices...."<<endl;

    std::string fNamedRdWReduced="dRdWReduced_"+np.str();
    adjIO_.writeMatrixBinary(dRdWReduced_,fNamedRdWReduced);
    if(debugMode) adjIO_.writeMatrixASCII(dRdWReduced_,fNamedRdWReduced);

    std::string fNamedRdFFDReduced="dRdFFDReduced_"+np.str();
    adjIO_.writeMatrixBinary(dRdFFDReduced_,fNamedRdFFDReduced);
    if(debugMode) adjIO_.writeMatrixASCII(dRdFFDReduced_,fNamedRdFFDReduced);

    Info<<"Writing the reduced matrices...."<<endl;
    Info<<"Done! "<<runTime_.elapsedClockTime()<<" s"<<endl;

}

void ReducedOrderModeling::perturbStatesMF(label n)
{
    // now we need to perturb eps*v to w, where v is the nth column of the phi matrix

    // first extract v: the nth column of the phiMat
    Vec colVec;
    VecCreate(PETSC_COMM_WORLD,&colVec);
    VecSetSizes(colVec,localSize_,PETSC_DECIDE);
    VecSetFromOptions(colVec);
    VecZeroEntries(colVec);
    MatGetColumnVector(svdPhiWMat_,colVec,n);

    PetscScalar *colVecArray;
    VecGetArray(colVec,&colVecArray);

    // perturb volVectorStates, modified based on perturbStates in adjDev
    const objectRegistry& db_(mesh_.thisDb());
    forAll(adjReg_.volVectorStates,idxI)                                           
    {        
        // create state and stateRef
        makeState(volVectorStates,volVectorField,adjReg_); 
        makeStateRef(volVectorStates,volVectorField,adjReg_);                 
                                                                               
        forAll(mesh_.cells(),cellI)                                            
        {                                                                      
            // check if this state's color = coloI
            for(label i=0;i<3;i++)                                         
            {
                PetscInt local = adjIdx_.getLocalAdjointStateIndex(stateName,cellI,i);
                scalar val = colVecArray[local];
                state[cellI][i] = stateRef[cellI][i] + mfStep*val;                    
            }                                                                  
        } 
        // correct BC
        state.correctBoundaryConditions();                                                                    
    }

    // perturb volScalarStates
    forAll(adjReg_.volScalarStates,idxI)
    {
        // create state and stateRef
        makeState(volScalarStates,volScalarField,adjReg_);
        makeStateRef(volScalarStates,volScalarField,adjReg_);
        
        forAll(mesh_.cells(),cellI)
        {         
            PetscInt local = adjIdx_.getLocalAdjointStateIndex(stateName,cellI);
            scalar val = colVecArray[local];
            state[cellI] = stateRef[cellI] + mfStep*val;    
        }
        
        // correct BC
        state.correctBoundaryConditions();  
    }

    // perturb turbStates
    forAll(adjRAS_.turbStates,idxI)
    {
        // create state and stateRef
        makeState(turbStates,volScalarField,adjRAS_);
        makeStateRef(turbStates,volScalarField,adjRAS_);
        
        forAll(mesh_.cells(),cellI)
        {         
            PetscInt local = adjIdx_.getLocalAdjointStateIndex(stateName,cellI);
            scalar val = colVecArray[local];
            state[cellI] = stateRef[cellI] + mfStep*val;    
        }

    }
    // BC for turbStates are implemented in the AdjRAS class
    adjRAS_.correctTurbBoundaryConditions();

    // perturb surfaceScalarStates
    forAll(adjReg_.surfaceScalarStates,idxI)
    {
        // create state and stateRef
        makeState(surfaceScalarStates,surfaceScalarField,adjReg_);
        makeStateRef(surfaceScalarStates,surfaceScalarField,adjReg_);
        
        forAll(mesh_.faces(),faceI)
        {
            // check if this state's color = coloI
            PetscInt local = adjIdx_.getLocalAdjointStateIndex(stateName,faceI);
            scalar val = colVecArray[local]; 
            if (faceI < adjIdx_.nLocalInternalFaces)
            {
                state[faceI] = stateRef[faceI] + mfStep*val;   
            }
            else
            {
                label relIdx=faceI-adjIdx_.nLocalInternalFaces;
                label patchIdx=adjIdx_.bFacePatchI[relIdx];
                label faceIdx=adjIdx_.bFaceFaceI[relIdx];
                state.boundaryFieldRef()[patchIdx][faceIdx] = 
                    stateRef.boundaryField()[patchIdx][faceIdx] + mfStep*val;   
            }           
        }
    }

    // NOTE: we also need to update states that are related to the adjoint states but not
    // perturbed here. For example, in buoyantBoussinesqFoam, p is related to p_rgh;
    // however, we only perturb p_rgh in this function. To calculate the perturbed
    // p due to the p_rgh perturbation for calculating force, we need to do p=p_rgh+rhok*gh
    // Similar treatment is needed for rhok and alphat. Basically, any variables apprear in flow residual
    // calculation or objection function calculation that are not state variables need to be updated.
    // This function is implemented in child class
    adjDev_.updateIntermediateVariables();


    VecRestoreArray(colVec,&colVecArray);

    return;
}

void ReducedOrderModeling::setdRdWPhiMat(Mat matIn,label n)
{
    PetscInt    Istart, Iend;
    // get the local ownership range
    MatGetOwnershipRange(matIn,&Istart,&Iend);

    for(PetscInt j=Istart; j<Iend; j++)
    {
        label localIdx = adjIdx_.globalAdjointStateNumbering.toLocal(j);
        PetscScalar val = adjDev_.adjStateLocalIdx2PartDerivVal(localIdx);
        MatSetValues(matIn,1,&j,1,&n,&val,INSERT_VALUES);
    }
}

void ReducedOrderModeling::calcdRdWPhiMF(Mat dRdWPhi)
{
    MatZeroEntries(dRdWPhi);
    // compute dRdW*phi
    // do matrix-free for dRdW*Phi = [ R(w+v*mfSetp) - R(w) ] / mfStep
    label isRef=1, isPC=0;
    adjDev_.copyStates("Ref2Var");
    adjDev_.calcResiduals(isRef,isPC);
    adjRAS_.calcTurbResiduals(isRef,isPC);
    nFuncEvals_++;
    for(label nn=0;nn<nSamples;nn++)
    {
        this->perturbStatesMF(nn);

        adjDev_.calcResPartDeriv(mfStep,isPC);
        nFuncEvals_++;

        // reset perturbation
        adjDev_.copyStates("Ref2Var");

        this->setdRdWPhiMat(dRdWPhi,nn);

    }

    MatAssemblyBegin(dRdWPhi,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dRdWPhi,MAT_FINAL_ASSEMBLY);
}

void ReducedOrderModeling::calcReducedMatsMF()
{

    // ********************** Ar *********************
    // first compute dRdW*phi
    Info<<"Computing dRdW*Phi..."<<endl;

    this->calcdRdWPhiMF(dRdWPhi_);

    adjDev_.calcFlowResidualStatistics("verify");

    if(debugMode) 
    {
        adjIO_.writeMatrixASCII(dRdWPhi_,"dRdWPhi");
        adjIO_.writeMatrixBinary(dRdWPhi_,"dRdWPhi");
    }

    // now compute phiT*dRdW*Phi
    //Mat phiT;
    //MatCreateTranspose(svdPhiWMat_,&phiT);
    void initializedRdWMatReduced();
    Info<< "Computing phiT*dRdW*phi" << endl;
    MatTransposeMatMult(svdPhiWMat_,dRdWPhi_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdWReduced_);

    // ********************** Br *********************
    Info<<"Calculating dRdFFD... "<<endl;
    this->initializedRdFFDMat();
    adjDev_.calcdRdFFD(dRdFFD_);
    Info<< "Computing phiT*dRdFFD" << endl;
    void initializedRdFFDMatReduced();
    MatTransposeMatMult(svdPhiWMat_,dRdFFD_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdFFDReduced_);

    MatDestroy(&dRdFFD_);
    
    return;

}

void ReducedOrderModeling::setNewField(Vec vecIn,word mode)
{
    const objectRegistry& db_(mesh_.thisDb());

    label Istart,Iend;
    label cellI,comp,faceI;
    VecGetOwnershipRange(vecIn,&Istart,&Iend);

    PetscScalar *vecInArray;
    VecGetArray(vecIn,&vecInArray);

    for(label i=Istart;i<Iend;i++)
    {
        label relIdx = i-Istart;
        word stateName = adjIdx_.adjStateName4LocalAdjIdx[relIdx];
        scalar cellIFaceI = adjIdx_.cellIFaceI4LocalAdjIdx[relIdx];
        const word& stateType = adjIdx_.adjStateType[stateName];

        if(stateType == "volVectorState")
        {
            volVectorField& state
            (
                const_cast<volVectorField&>                                             
                ( 
                    db_.lookupObject<volVectorField>(stateName)
                )
            );

            cellI = round(cellIFaceI);
            comp = round(10*(cellIFaceI-cellI));
            if (mode=="add") state[cellI][comp] += vecInArray[relIdx];
            else if (mode=="set") state[cellI][comp] = vecInArray[relIdx];
            else FatalErrorIn("")<<"stateType not known"<< abort(FatalError);
            
        }
        else if(stateType == "volScalarState" or stateType == "turbState")
        {
            volScalarField& state
            (
                const_cast<volScalarField&>                                             
                ( 
                    db_.lookupObject<volScalarField>(stateName)
                )
            );

            cellI = round(cellIFaceI);
            if (mode=="add") state[cellI] += vecInArray[relIdx];
            else if (mode=="set") state[cellI] = vecInArray[relIdx];
            else FatalErrorIn("")<<"stateType not known"<< abort(FatalError);
        }
        else if(stateType == "surfaceScalarState")
        {
            surfaceScalarField& state
            (
                const_cast<surfaceScalarField&>                                             
                ( 
                    db_.lookupObject<surfaceScalarField>(stateName)
                )
            );
            faceI = round(cellIFaceI);
            if(faceI<adjIdx_.nLocalInternalFaces)
            {
                if (mode=="add") state[faceI] += vecInArray[relIdx];
                else if (mode=="set") state[faceI] = vecInArray[relIdx];
                else FatalErrorIn("")<<"stateType not known"<< abort(FatalError);
            }
            else
            {
                label relIdx=faceI-adjIdx_.nLocalInternalFaces;
                label patchIdx=adjIdx_.bFacePatchI[relIdx];
                label faceIdx=adjIdx_.bFaceFaceI[relIdx];
                if (mode=="add") state.boundaryFieldRef()[patchIdx][faceIdx] += vecInArray[relIdx];
                else if (mode=="set") state.boundaryFieldRef()[patchIdx][faceIdx] = vecInArray[relIdx];
                else FatalErrorIn("")<<"stateType not known"<< abort(FatalError);
            }
        }
        else
        {
            FatalErrorIn("")<<"stateType not known"<< abort(FatalError);
        }
    }
    
    VecRestoreArray(vecIn,&vecInArray);

    // need to update BCs and intermediate variable
    adjDev_.updateStateVariableBCs();
    adjDev_.updateIntermediateVariables();
}

void ReducedOrderModeling::writeNewField(word postfix)
{
    const objectRegistry& db_(mesh_.thisDb());
    // write
    forAll(adjReg_.volVectorStates,idxI)                                           
    {        
        // create state and stateRef
        makeState(volVectorStates,volVectorField,adjReg_); 
        word oldName=state.name();
        word newName=state.name()+postfix;
        state.correctBoundaryConditions();
        state.rename(newName);
        state.write();       
        state.rename(oldName);                                                             
    }

    forAll(adjReg_.volScalarStates,idxI)
    {
        // create state and stateRef
        makeState(volScalarStates,volScalarField,adjReg_);
        word oldName=state.name();
        word newName=state.name()+postfix;
        state.correctBoundaryConditions();
        state.rename(newName);
        state.write();       
        state.rename(oldName);   
    }

    forAll(adjRAS_.turbStates,idxI)
    {
        // create state and stateRef
        makeState(turbStates,volScalarField,adjRAS_);
        word oldName=state.name();
        word newName=state.name()+postfix;
        state.correctBoundaryConditions();
        state.rename(newName);
        state.write();       
        state.rename(oldName); 

    }

    // write nutROM 
    volScalarField& nut 
    (
        const_cast<volScalarField&>
        (
            db_.lookupObject<volScalarField>("nut")
        )
    );
    nut.rename("nut"+postfix);
    nut.write();
    nut.rename("nut");

    forAll(adjReg_.surfaceScalarStates,idxI)
    {
        // create state and stateRef
        makeState(surfaceScalarStates,surfaceScalarField,adjReg_);
        word oldName=state.name();
        word newName=state.name()+postfix;
        state.rename(newName);
        state.write();       
        state.rename(oldName); 
    }

}

void ReducedOrderModeling::solveOnlineLinear()
{
    this->initializeOnlineLinear();

    KSP ksp;
    // add options and initialize ksp
    dictionary adjOptions;
    adjOptions.add("GMRESRestart",adjIO_.nkGMRESRestart);
    adjOptions.add("GlobalPCIters",adjIO_.nkGlobalPCIters);
    adjOptions.add("ASMOverlap",adjIO_.nkASMOverlap);
    adjOptions.add("LocalPCIters",adjIO_.nkLocalPCIters);
    adjOptions.add("JacMatReOrdering",adjIO_.nkJacMatReOrdering);
    adjOptions.add("PCFillLevel",adjIO_.nkPCFillLevel);
    adjOptions.add("GMRESMaxIters",adjIO_.nkGMRESMaxIters);
    adjOptions.add("GMRESRelTol",adjIO_.adjGMRESRelTol);
    adjOptions.add("GMRESAbsTol",adjIO_.adjGMRESAbsTol);
    adjOptions.add("printInfo",1);
    //Info<<adjOptions<<endl;
    adjDev_.createMLRKSP(&ksp,dRdWReduced_,dRdWReduced_,adjOptions);

    Vec RHS,deltaWVecReduced, deltaWVec;

    VecCreate(PETSC_COMM_WORLD,&RHS);
    VecSetSizes(RHS,PETSC_DECIDE,nSamples);
    VecSetFromOptions(RHS);
    VecDuplicate(RHS,&deltaWVecReduced);

    VecZeroEntries(deltaWVecReduced);
    VecZeroEntries(RHS);
 
    VecCreate(PETSC_COMM_WORLD,&deltaWVec);
    VecSetSizes(deltaWVec,localSize_,PETSC_DECIDE);
    VecSetFromOptions(deltaWVec);
    VecZeroEntries(deltaWVec);

    MatMult(dRdFFDReduced_,deltaFFDVec_,RHS);
    VecScale(RHS,-1.0);

    KSPSolve(ksp,RHS, deltaWVecReduced);
    //Print convergence information
    label its;
    KSPGetIterationNumber(ksp,&its);
    PetscScalar finalResNorm;
    KSPGetResidualNorm(ksp,&finalResNorm);
    PetscPrintf
    (
        PETSC_COMM_WORLD,
        "Main iteration %D KSP Residual norm %14.12e %d s \n",
        its,
        finalResNorm,
        runTime_.elapsedClockTime()
    );
    PetscPrintf(PETSC_COMM_WORLD,"Total iterations %D\n",its);

    MatMult(svdPhiWMat_,deltaWVecReduced,deltaWVec);

    if(debugMode)
    {
        adjIO_.writeVectorASCII(deltaWVecReduced,"deltaWVecReduced");
        adjIO_.writeVectorASCII(deltaWVec,"deltaWVec");
    }

    this->setNewField(deltaWVec,"add");
    adjObj_.writeObjFuncValues();
    this->writeNewField("ROM");

    // now write the modes based on the phiMat_
    // extract the nth column of the phiMat
    Vec colVec;
    VecCreate(PETSC_COMM_WORLD,&colVec);
    VecSetSizes(colVec,localSize_,PETSC_DECIDE);
    VecSetFromOptions(colVec);

    for(label n=0;n<nSamples;n++)
    {
        VecZeroEntries(colVec);
        MatGetColumnVector(svdPhiWMat_,colVec,n);
        this->setNewField(colVec,"set");
        word fieldPostfix="Mode"+Foam::name(n);
        this->writeNewField(fieldPostfix);
    }

    Info<<"Done! "<<runTime_.elapsedClockTime()<<" s"<<endl;

}

void ReducedOrderModeling::solveOnlineNonlinear()
{
    this->initializeOnlineNonlinear();

    // we use the Newton method to solver rResiduals(rStates)=0
    this->solveNK();

    Info<<"Done! "<<runTime_.elapsedClockTime()<<" s"<<endl;

}

void ReducedOrderModeling::solveNK()
{
    Info<<endl<<"***** Solving ROM Nonlinear Equation Using Newton-Krylov Method *****"<<endl;

    // KSP
    KSP ksp;
    // number of GMRES linear iterations
    label GMRESIters;
    // EW parameters
    scalar rVecNorm;
    scalar rTol=romNKGMRESRTol;
    scalar aTol=1.0e-12;
    // reduced Jacobians
    Mat rdRdW,rdRdWPC;
    // total residual norm including all variables
    scalar totalResNorm0;
    scalar totalResNorm;
    // these vars are for store the tolerance for GMRES linear solution
    scalar rGMRESHist[adjIO_.nkGMRESMaxIters+1];
    label nGMRESIters=adjIO_.nkGMRESMaxIters+1;
    // these vars are for store the FD step size for GMRES linear solution
    scalar nkHHist[adjIO_.nkGMRESMaxIters+1];
    label nkHN=adjIO_.nkGMRESMaxIters+1;
    for(label i=0;i<adjIO_.nkGMRESMaxIters+1;i++) nkHHist[i]=0.0;

    // create reduced Jacobian and set function evaluation 
    MatCreateMFFD(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,nSamples,nSamples,&rdRdW);
    MatMFFDSetFunction(rdRdW,FormFunction,this);

    if (romNKMFFDH > 0.0) MatMFFDSetCheckh(rdRdW,ComputeMFFDH,this);

    // initialize and calculate rdRdWPC
    this->initializeReducedJacobian(&rdRdWPC);
    // save the unperturb residual statistics as reference
    // we will verify against this reference to ensure consistent residuals after perturbing
    // and resetting states
    adjDev_.calcFlowResidualStatistics("set");
    adjDev_.calcFlowResidualStatistics("print");
    //this->calcReducedJacobian(rdRdWPC); 
    adjIO_.readMatrixBinary(rdRdWPC,"rdRdWPC");
    if (debugMode) adjIO_.writeMatrixASCII(rdRdWPC,"rdRdWPC");

    // first compute/assign the initial wVec and rVec and compute the initial totalResNorm0
    VecZeroEntries(wVecReduced_);
    VecZeroEntries(rVecReduced_);  
    VecZeroEntries(wVecFull_);  
    this->NKSetVecs(wVecFull_,"Var2Vec",1.0,"");
    MatMultTranspose(svdPhiWMat_,wVecFull_,wVecReduced_);
    if (useLSPG) this->calcdRdWPhiMF(dRdWPhi_); // we need to compute the initial dRdWPhi for LSPG
    this->NKCalcResidualsReduced(wVecReduced_,rVecReduced_);

    // compute the initial norm
    VecNorm(rVecReduced_,NORM_2,&totalResNorm0);
    totalResNorm=totalResNorm0;
    scalar totalResNormFull=this->getResNorm("total");
    scalar turbResNormFull=this->getResNorm("turb");

    // add options and initialize ksp
    dictionary adjOptions;
    adjOptions.add("GMRESRestart",100);
    adjOptions.add("GlobalPCIters",0);
    adjOptions.add("ASMOverlap",1);
    adjOptions.add("LocalPCIters",1);
    adjOptions.add("JacMatReOrdering","rcm");
    adjOptions.add("PCFillLevel",0);
    adjOptions.add("GMRESMaxIters",100);
    adjOptions.add("GMRESRelTol",rTol);
    adjOptions.add("GMRESAbsTol",aTol);
    adjOptions.add("printInfo",0);
    //Info<<adjOptions<<endl;
    adjDev_.createMLRKSP(&ksp,rdRdW,rdRdWPC,adjOptions);

    // initialize objFuncs
    HashTable<scalar> objFuncs;
    forAll(adjIO_.objFuncs,idxI)
    {
        word objFunc = adjIO_.objFuncs[idxI];
        adjObj_.calcObjFuncs(objFunc,0);
        scalar val = adjObj_.getObjFunc(objFunc);
        objFuncs.set(objFunc,val);
    }

    // print the initial convergence information
    this->printConvergenceInfo("printHeader",objFuncs);
    this->printConvergenceInfo
    (
        "printConvergence",
        objFuncs,
        0,
        0,
        "  NK  ",
        1.0,
        0.0,
        1.0,
        scalar(adjDev_.getRunTime()),
        turbResNormFull,
        totalResNormFull,
        totalResNorm
    );

    // main loop for NK
    Vec rVecReducedBase,rhs, dWVecReduced,rVecReducedNew,wVecReducedNew;
    VecDuplicate(wVecReduced_,&dWVecReduced);
    VecDuplicate(rVecReduced_,&rhs);
    VecDuplicate(rVecReduced_,&rVecReducedBase);
    VecDuplicate(rVecReduced_,&rVecReducedNew);
    VecDuplicate(wVecReduced_,&wVecReducedNew);
    for(label iterI=1;iterI<romNKMaxIts;iterI++)
    {
        // check if the presribed tolerances are met
        if (totalResNorm<romNKAbsTol)
        {
            Info<<"Absolute Tolerance "<<totalResNorm<<" less than the presribed nkAbsTol "<<romNKAbsTol<<endl;
            Info<<"NK completed!"<<endl;
            break;
        }

        // set up rGMRESHist to save the tolerance history for the GMRES solution
        KSPSetResidualHistory(ksp,rGMRESHist,nGMRESIters,PETSC_TRUE);
        // set up the FD step history and print it for debugging
        MatMFFDSetHHistory(rdRdW,nkHHist,nkHN);

        // before solving the ksp, form the baseVector for matrix-vector products.
        // Note: we need to apply normalize-states to the baseVector
        // Note that we also scale the dRdW*psi 
        // in AdjointNewtonKrylov::FormFunction
        VecCopy(rVecReduced_,rVecReducedBase); 
        //this->setNormalizeStatesScaling2Vec(rVecBase);
        MatMFFDSetBase(rdRdW,wVecReduced_,rVecReducedBase);

        // for each major iteration, update the dRdWPhi matrix for LSPG
        if (useLSPG) this->calcdRdWPhiMF(dRdWPhi_);

        // solve the linear system
        // we should use rVec0 as the rhs, however, we need to normalize
        // the states, so we create this temporary rhs vec to store 
        // the scaled rVec. Note that we also scale the dRdW*psi 
        // in AdjointNewtonKrylov::FormFunction
        VecCopy(rVecReduced_,rhs);
        KSPSolve(ksp,rhs,dWVecReduced);

        // get linear solution rTol: linRes
        KSPGetIterationNumber(ksp,&GMRESIters);
        scalar linRes=rGMRESHist[GMRESIters]/rGMRESHist[0];

        // do a line search and update states
        VecZeroEntries(rVecReducedNew);
        VecZeroEntries(wVecReducedNew);
        scalar stepSize=this->NKLineSearch(wVecReduced_,rVecReduced_,dWVecReduced,wVecReducedNew,rVecReducedNew);
        VecCopy(rVecReducedNew,rVecReduced_);
        VecCopy(wVecReducedNew,wVecReduced_);

        // update the norm for printting convergence info
        VecNorm(rVecReduced_,NORM_2,&rVecNorm);
        totalResNorm=rVecNorm;
        scalar totalResNormFull=this->getResNorm("total");
        scalar turbResNormFull=this->getResNorm("turb");
    
        // update objective function values
        forAll(adjIO_.objFuncs,idxI)
        {
            word objFunc = adjIO_.objFuncs[idxI];
            adjObj_.calcObjFuncs(objFunc,0);
            scalar val = adjObj_.getObjFunc(objFunc);
            objFuncs.set(objFunc,val);
        }

        // print convergence info
        this->printConvergenceInfo
        (
            "printConvergence",
            objFuncs,
            iterI,
            nFuncEvals_,
            "  NK  ",
            stepSize,
            linRes,
            1.0,
            scalar(adjDev_.getRunTime()),
            turbResNormFull,
            totalResNormFull,
            totalResNorm
        );

        if (debugMode)
        {
            Info<<"NK Iter: "<<iterI<<endl;
            for(label i=0;i<nkHN;i++)
            {
                if(fabs(nkHHist[i])>1e-16)
                {
                    Info<<"GMRES iter"<<i<<" FD step: "<<nkHHist[i]<<endl;
                }
            }
            std::ostringstream nn("");
            nn<<iterI;
            std::string name1="rVecReduced_"+nn.str();
            std::string name2="wVecReduced_"+nn.str();
            std::string name3="rVecFull_"+nn.str();
            std::string name4="wVecFull_"+nn.str();
            adjIO_.writeVectorASCII(rVecReduced_,name1);
            adjIO_.writeVectorASCII(wVecReduced_,name2);
            adjIO_.writeVectorASCII(rVecFull_,name3);
            adjIO_.writeVectorASCII(wVecFull_,name4);
        }

    }

    // assign the latest wVec to variables
    VecZeroEntries(wVecFull_);
    MatMult(svdPhiWMat_,wVecReduced_,wVecFull_);
    this->NKSetVecs(wVecFull_,"Vec2Var",1.0,"");
    adjObj_.writeObjFuncValues();
    this->writeNewField("ROM");

    adjDev_.calcFlowResidualStatistics("print");
    
}


scalar ReducedOrderModeling::NKLineSearch
(
    const Vec wVec,
    const Vec rVec,
    const Vec dWVec,
    Vec wVecNew,
    Vec rVecNew
)
{
    // basic non monotonic line search using backtracking 
    
    scalar gamma = 0.999; // line search sufficient decrease coefficient
    scalar sigma = 0.5;    // line search reduce factor
    scalar alpha = 1.0;      // initial step size

    PetscErrorCode ierr;

    // compute the norms
    scalar rVecNorm=0.0, rVecNewNorm=0.0;
    VecNorm(rVec,NORM_2,&rVecNorm);

    // initial step
    alpha=1.0;

    // actual backtracking
    for(label i=0;i<romNKGMRESMaxLS;i++)
    {
        // compute new w value wNew = w-dW Note: dW is the solution from KSP it is dW=wOld-wNew
        VecWAXPY(wVecNew,-alpha,dWVec,wVec);

        // using the current state to compute rVecNew
        this->NKCalcResidualsReduced(wVecNew,rVecNew);

        // compute the rVecNorm at new w
        ierr=VecNorm(rVecNew,NORM_2,&rVecNewNorm);

        if(ierr == PETSC_ERR_FP)
        {
            // floating point seg fault, just reduce the step size
            alpha = alpha*sigma;
            Info<<"Seg Fault in line search"<<endl;
            continue;
        }
        else if(this->checkNegativeTurb())
        {
            // have negative turbulence, reduce the step size
            alpha = alpha*sigma;
            Info<<"Negative turbulence in line search"<<endl;
            continue;
        }
        else if (rVecNewNorm <= rVecNorm*gamma)  // Sufficient reduction
        {
            return alpha;
        }
        else // reduce step
        {
            alpha = alpha * sigma;
        }

    }

    return alpha;

}

label ReducedOrderModeling::checkNegativeTurb()
{
    if(adjIO_.nkSegregatedTurb) return 0;

    const objectRegistry& db=mesh_.thisDb();
    forAll(adjIdx_.adjStateNames,idxI)
    {
        word stateName = adjIdx_.adjStateNames[idxI];
        word stateType = adjIdx_.adjStateType[stateName];
        if (stateType == "turbState")
        {
            const volScalarField& var =  db.lookupObject<volScalarField>(stateName) ;
            forAll(var,idxJ)
            {
                if (var[idxJ]<0)
                {
                    return 1;
                }
            }
            forAll(var.boundaryField(),patchI)
            {
                forAll(var.boundaryField()[patchI],faceI)
                {
                    if (var.boundaryField()[patchI][faceI]<0)
                    {
                        return 1;
                    }
                }
            }
        }
    }

    return 0;
}


scalar ReducedOrderModeling::getEWTol
(
    scalar norm, 
    scalar oldNorm, 
    scalar rTolLast
)
{
    // There are the default EW Parameters from PETSc. They seem to work well
    // version:  2
    // rTolLast: 0.1
    // rTolMax:  0.9
    // gamma:    1.0
    // alpha:    1.61803398874989
    // threshold: 0.1

    scalar rTolMax   = adjIO_.nkEWRTolMax;
    scalar gamma     = 1.0;
    scalar alpha     = (1.0+Foam::sqrt(5.0))/2.0;
    scalar threshold = 0.10;
    // We use version 2:
    scalar rTol = gamma*Foam::pow(norm/oldNorm,alpha);
    scalar sTol = gamma*Foam::pow(rTolLast,alpha);

    if (sTol > threshold)
    {
       rTol = max(rTol, sTol);
    }

    // Safeguard: avoid rtol greater than one
    rTol = min(rTol, rTolMax);

    return rTol;
}

void ReducedOrderModeling::calcReducedJacobian(Mat matIn)
{

    Info<<"Computing reduced Jacobian "<<endl;
    
    PetscInt Istart, Iend;
    MatGetOwnershipRange(matIn,&Istart,&Iend);

    // comptue wVecFullRef and rVecFullRef
    // assign the current reference states from OpenFOAM to wVecFull_
    Vec wVecFullRef;
    VecDuplicate(wVecFull_,&wVecFullRef);
    VecZeroEntries(wVecFullRef);
    this->NKSetVecs(wVecFullRef,"Var2Vec",1.0,"");
    // compute the referene full residual rVecFull
    Vec rVecFullRef;
    VecDuplicate(rVecFull_,&rVecFullRef);
    VecZeroEntries(rVecFullRef);
    this->NKCalcResidualsFull(wVecFullRef,rVecFullRef);
    // then compute rVecReducedRef
    Vec rVecReducedRef;
    VecDuplicate(rVecReduced_,&rVecReducedRef);
    VecZeroEntries(rVecReducedRef);
    if(useLSPG)
    {
        this->calcdRdWPhiMF(dRdWPhi_);
        MatMultTranspose(dRdWPhi_,rVecFullRef,rVecReducedRef);
    }
    else
    {
        //MatMultTranspose(svdPhiRMat_,rVecFullRef,rVecReducedRef);
        MatMultTranspose(svdPhiWMat_,rVecFullRef,rVecReducedRef);
    }
    
    Vec wVecReducedDelta, wVecFullDelta;
    VecDuplicate(wVecReduced_,&wVecReducedDelta);
    VecDuplicate(wVecFull_,&wVecFullDelta);

    // get the reference wVecFull_
    //this->NKSetVecs(wVecFull_,"Var2Vec",1.0,"");

    // perturb
    for(PetscInt nn=0;nn<nSamples;nn++)
    {
        // compute perturbed reduced wVec
        VecZeroEntries(wVecReducedDelta);
        scalar deltaVal=0.001;
        VecSetValue(wVecReducedDelta,nn,deltaVal,INSERT_VALUES);
        VecAssemblyBegin(wVecReducedDelta);
        VecAssemblyEnd(wVecReducedDelta);
        // compute perturbed full wVec
        VecZeroEntries(wVecFullDelta);
        MatMult(svdPhiWMat_,wVecReducedDelta,wVecFullDelta);

        // now the new wVecFull=wVecFull+wVecFullDelta
        VecZeroEntries(wVecFull_);
        VecWAXPY(wVecFull_,1.0,wVecFullRef,wVecFullDelta);

        // compute the perturbed full residual
        VecZeroEntries(rVecFull_);
        this->NKCalcResidualsFull(wVecFull_,rVecFull_);

        // then compute rVecReduced_
        VecZeroEntries(rVecReduced_);
        if(useLSPG)
        {
            this->calcdRdWPhiMF(dRdWPhi_);
            MatMultTranspose(dRdWPhi_,rVecFull_,rVecReduced_);
        }
        else
        {
            //MatMultTranspose(svdPhiRMat_,rVecFull_,rVecReduced_);
            MatMultTranspose(svdPhiWMat_,rVecFull_,rVecReduced_);
        }
        
        // now we know rVecReducedRef and rVecReduced_, we can use FD to compute partials
        VecAXPY(rVecReduced_,-1.0,rVecReducedRef);
        VecScale(rVecReduced_,1.0/deltaVal);
        // assign it to matIn
        scalar vecVal;
        for (PetscInt idx=Istart;idx<Iend;idx++)
        {
            VecGetValues(rVecReduced_,1,&idx,&vecVal);
            MatSetValues(matIn,1,&idx,1,&nn,&vecVal,INSERT_VALUES);
        }

    }
    MatAssemblyBegin(matIn,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matIn,MAT_FINAL_ASSEMBLY);
}

void ReducedOrderModeling::initializeReducedJacobian(Mat* PCMat)
{
   
    // initialize PCMat
    MatCreate(PETSC_COMM_WORLD,PCMat);
    MatSetSizes(*PCMat,PETSC_DECIDE,PETSC_DECIDE,nSamples,nSamples);
    MatSetFromOptions(*PCMat);
    MatMPIAIJSetPreallocation(*PCMat,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(*PCMat,nSamples,NULL);
    MatSetUp(*PCMat);
    MatZeroEntries(*PCMat);
    

}

PetscErrorCode ReducedOrderModeling::FormFunction(void* ctx,Vec wVec,Vec rVec)
{
    
    // Given an input vec wVec, calculate the function vec rWec
    // wVec: vector storing all the states
    // rVec: vector storing the residuals

    ReducedOrderModeling *rom = (ReducedOrderModeling*) ctx;

    rom->NKCalcResidualsReduced(wVec,rVec);

    //rom->setNormalizeStatesScaling2Vec(rVec);
    
    return 0;
    
}

PetscErrorCode ReducedOrderModeling::ComputeMFFDH(void* ctx,Vec vec1,Vec vec2, PetscScalar* h)
{
    
    
    ReducedOrderModeling *rom = (ReducedOrderModeling*) ctx;

    *h=rom->romNKMFFDH;
    
    return 0;
    
}

void ReducedOrderModeling::NKCalcResidualsReduced(Vec wVec,Vec rVec)
{
    // Given an input vec wVec, calculate the function vec rWec
    // wVec: vector storing all the states
    // rVec: vector storing the residuals

    // NOTE: the reduced residual is   Phi^T * R = 0
    
    // first compute the full length state vector from the reduced state vector
    VecZeroEntries(wVecFull_);
    MatMult(svdPhiWMat_,wVec,wVecFull_);

    // assign the wVecFull to the OpenFOAM field variables
    this->NKSetVecs(wVecFull_,"Vec2Var",1.0,"");

    // compute the full length residual in DAFoam
    label isRef=0,isPC=0;
    adjDev_.updateStateVariableBCs();
    adjDev_.calcResiduals(isRef,isPC);
    adjRAS_.calcTurbResiduals(isRef,isPC);
    nFuncEvals_++;
    //Info<<"Calling NKCalcResidualsReduced"<<endl;

    // assign the OpenFOAM field variables back to vVecFull_
    VecZeroEntries(rVecFull_);
    this->NKSetVecs(rVecFull_,"Var2Vec",1.0,"Res");

    // convert the full length residual to the reduced residual vector
    VecZeroEntries(rVec);
    if(useLSPG)
    {
        MatMultTranspose(dRdWPhi_,rVecFull_,rVec);
    }
    else
    {
        //MatMultTranspose(svdPhiRMat_,rVecFull_,rVec);
        MatMultTranspose(svdPhiWMat_,rVecFull_,rVec);
    }
    
    return;

}

void ReducedOrderModeling::NKCalcResidualsFull(Vec wVec,Vec rVec)
{
    // Given an input vec wVec, calculate the function vec rWec
    // wVec: vector storing all the states
    // rVec: vector storing the residuals

    // assign the wVecFull to the OpenFOAM field variables
    this->NKSetVecs(wVec,"Vec2Var",1.0,"");

    // compute the full length residual in DAFoam
    label isRef=0,isPC=0;
    adjDev_.updateStateVariableBCs();
    adjDev_.calcResiduals(isRef,isPC);
    adjRAS_.calcTurbResiduals(isRef,isPC);
    nFuncEvals_++;
    //Info<<"Calling NKCalcResidualsFull"<<endl;

    // assign the OpenFOAM field variables back to vVecFull_
    VecZeroEntries(rVec);
    this->NKSetVecs(rVec,"Var2Vec",1.0,"Res");

    return;

}

void ReducedOrderModeling::NKSetVecs
(
    Vec vecX, 
    word mode, 
    scalar scaleFactor, 
    word postFix
)
{
    
    // this is a general function to assign vec to/from states/residuals

    PetscInt    Istart, Iend;
    VecGetOwnershipRange(vecX,&Istart,&Iend);
    PetscScalar* vecXArray;
    const PetscScalar* vecXArrayConst;
    if(mode == "Vec2Var" || mode == "VecAdd2Var")
    {
        VecGetArrayRead(vecX,&vecXArrayConst);
    }
    else if (mode == "Var2Vec" || mode == "VarAdd2Vec")
    {
        VecGetArray(vecX,&vecXArray);
    }
    else
    {
        FatalErrorIn("")<<"mode not valid"<< abort(FatalError);
    }
    
    for(PetscInt i=Istart;i<Iend;i++)
    {
        label localIdx =  i-Istart;

        word stateName = adjIdx_.adjStateName4LocalAdjIdx[localIdx];
        word varName = stateName+postFix;
        const word& stateType = adjIdx_.adjStateType[stateName];
        
        scalar cellIFaceI =  adjIdx_.cellIFaceI4LocalAdjIdx[localIdx];

        if(stateType == "volVectorState")
        {

            volVectorField& var = const_cast<volVectorField&>( db_.lookupObject<volVectorField>(varName) );
            label cellI,comp;
            cellI = round(cellIFaceI);
            comp = round(10*(cellIFaceI-cellI));

            if(mode == "Vec2Var")
            {
                var[cellI][comp] = scaleFactor*vecXArrayConst[localIdx];
            }
            else if (mode == "Var2Vec")
            {
                vecXArray[localIdx] = scaleFactor*var[cellI][comp];
            }
            else if (mode == "VecAdd2Var")
            {
                var[cellI][comp] += scaleFactor*vecXArrayConst[localIdx];
            }
            else if (mode == "VarAdd2Vec")
            {
                vecXArray[localIdx] += scaleFactor*var[cellI][comp];
            }
            else
            {
                 FatalErrorIn("")<<"mode not valid"<< abort(FatalError);
            }
            
        }
        else if(stateType == "volScalarState")
        {

            volScalarField& var = const_cast<volScalarField&>( db_.lookupObject<volScalarField>(varName) );
            label cellI;
            cellI = round(cellIFaceI);

            if(mode == "Vec2Var")
            {
                var[cellI] = scaleFactor*vecXArrayConst[localIdx];
            }
            else if (mode == "Var2Vec")
            {
                vecXArray[localIdx] = scaleFactor*var[cellI];
            }
            else if (mode == "VecAdd2Var")
            {
                var[cellI] += scaleFactor*vecXArrayConst[localIdx];
            }
            else if (mode == "VarAdd2Vec")
            {
                vecXArray[localIdx] += scaleFactor*var[cellI];
            }
            else
            {
                 FatalErrorIn("")<<"mode not valid"<< abort(FatalError);
            }

        }
        else if(stateType == "turbState" and !adjIO_.nkSegregatedTurb)
        {

            volScalarField& var = const_cast<volScalarField&>( db_.lookupObject<volScalarField>(varName) );
            label cellI;
            cellI = round(cellIFaceI);

            if(mode == "Vec2Var")
            {
                var[cellI] = scaleFactor*vecXArrayConst[localIdx];
            }
            else if (mode == "Var2Vec")
            {
                vecXArray[localIdx] = scaleFactor*var[cellI];
            }
            else if (mode == "VecAdd2Var")
            {
                var[cellI] += scaleFactor*vecXArrayConst[localIdx];
            }
            else if (mode == "VarAdd2Vec")
            {
                vecXArray[localIdx] += scaleFactor*var[cellI];
            }
            else
            {
                 FatalErrorIn("")<<"mode not valid"<< abort(FatalError);
            }
        }
        else if (stateType == "surfaceScalarState")
        {

            surfaceScalarField& var = const_cast<surfaceScalarField&>( db_.lookupObject<surfaceScalarField>(varName) );
            label faceI=round(cellIFaceI);
            if(faceI<mesh_.nInternalFaces())
            {
                if(mode == "Vec2Var")
                {
                    var[faceI] = scaleFactor*vecXArrayConst[localIdx];
                }
                else if (mode == "Var2Vec")
                {
                    vecXArray[localIdx] = scaleFactor*var[faceI];
                }
                else if (mode == "VecAdd2Var")
                {
                    var[faceI] += scaleFactor*vecXArrayConst[localIdx];
                }
                else if (mode == "VarAdd2Vec")
                {
                    vecXArray[localIdx] += scaleFactor*var[faceI];
                }
                else
                {
                     FatalErrorIn("")<<"mode not valid"<< abort(FatalError);
                }
            }
            else
            {
                label relIdx=faceI-mesh_.nInternalFaces();
                label patchIdx=adjIdx_.bFacePatchI[relIdx];
                label faceIdx=adjIdx_.bFaceFaceI[relIdx];

                if(mode == "Vec2Var")
                {
                    var.boundaryFieldRef()[patchIdx][faceIdx] = scaleFactor*vecXArrayConst[localIdx];
                }
                else if (mode == "Var2Vec")
                {
                    vecXArray[localIdx] = scaleFactor*var.boundaryFieldRef()[patchIdx][faceIdx];
                }
                else if (mode == "VecAdd2Var")
                {
                    var.boundaryFieldRef()[patchIdx][faceIdx] += scaleFactor*vecXArrayConst[localIdx];
                }
                else if (mode == "VarAdd2Vec")
                {
                    vecXArray[localIdx] += scaleFactor*var.boundaryFieldRef()[patchIdx][faceIdx];
                }
                else
                {
                     FatalErrorIn("")<<"mode not valid"<< abort(FatalError);
                }

            }
        }
        else
        {
            FatalErrorIn("")<<"statetype not valid"<< abort(FatalError);
        }
    }
    
    if(mode == "Vec2Var" || mode == "VecAdd2Var")
    {
        VecRestoreArrayRead(vecX,&vecXArrayConst);
    }
    else if (mode == "Var2Vec" || mode == "VarAdd2Vec")
    {
        VecRestoreArray(vecX,&vecXArray);
        //Info<<"Pass NKFormFunctoin5"<<endl;
        VecAssemblyBegin(vecX);
        VecAssemblyEnd(vecX);
    }
    else
    {
        FatalErrorIn("")<<"mode not valid"<< abort(FatalError);
    }

    if(mode == "Vec2Var" or mode == "VecAdd2Var")
    {
        adjDev_.copyStates("Var2Ref"); // need to update the reference states as well
        adjDev_.updateStateVariableBCs();
    }
    
    return;

}


void ReducedOrderModeling::printConvergenceInfo
(
    word mode,
    HashTable<scalar> objFuncs,
    label mainIter,
    label nFuncEval,
    word solverType,
    scalar step,
    scalar linRes,
    scalar CFL,
    scalar runTime,
    scalar turbNorm,
    scalar phiNorm,
    scalar totalNorm
    
)
{
    if (mode=="printHeader")
    {
        word printInfo="";
        printInfo +=   "--------------------------------------------------------------------------------------------------------";
        forAll(objFuncs.toc(),idxI)
        {
            printInfo += "------------------------";
        }
        printInfo += "\n";

        printInfo +=   "| Main | Func | Solv |  Step  |   Lin   |   CFL   |   RunTime   |  Res Norm  |  Res Norm  |  Res Norm  |";
        forAll(objFuncs.toc(),idxI)
        {
            word key = objFuncs.toc()[idxI];
            label keySize=key.size();
            for(label i=0;i<13-keySize;i++) printInfo += " ";
            printInfo += key+"          |";
        }
        printInfo += "\n";

        printInfo +=   "| Iter | Eval | Type |        |   Res   |         |   (s)       |  (Turb)    |  (Full)    |  (Reduced) |";
        forAll(objFuncs.toc(),idxI)
        {
            word key = objFuncs.toc()[idxI];
            printInfo += "                       |";
        }
        printInfo += "\n";

        printInfo +=   "--------------------------------------------------------------------------------------------------------";
        forAll(objFuncs.toc(),idxI)
        {
            printInfo += "------------------------";
        }
        printInfo += "\n";
    
        PetscPrintf
        (
            PETSC_COMM_WORLD,
            printInfo.c_str()
        );
    }
    else if (mode=="printConvergence")
    {
        PetscPrintf
        (
            PETSC_COMM_WORLD,
            " %6d %6d %s  %.4f   %.1e   %.2e   %.4e   %.4e   %.4e   %.4e",
            mainIter,
            nFuncEval,
            solverType.c_str(),
            step,
            linRes,
            CFL,
            runTime,
            turbNorm,
            phiNorm,
            totalNorm
        );
        forAll(objFuncs.toc(),idxI)
        {
            word key = objFuncs.toc()[idxI];
            scalar val = objFuncs[key];
            PetscPrintf
            (
                PETSC_COMM_WORLD,
                "   %.15e",
                val
            );
        }
        PetscPrintf
        (
            PETSC_COMM_WORLD,
            "\n"
        );
    }
    else
    {
        FatalErrorIn("")<< "mode not found!"<<abort(FatalError);
    }
}

scalar ReducedOrderModeling::getResNorm(word mode)
{

    scalar totalResNorm2=0.0;
    scalar turbResNorm2=0.0;
    scalar phiResNorm2=0.0;

    forAll(adjReg_.volVectorStates,idxI)
    {
        const word stateName = adjReg_.volVectorStates[idxI];
        const word resName = stateName+"Res";                   
        const volVectorField& stateRes = db_.lookupObject<volVectorField>(resName); 
        
        vector vecResNorm2(0,0,0);
        forAll(stateRes,cellI)
        {
            vecResNorm2.x()+=Foam::pow(stateRes[cellI].x(),2.0);
            vecResNorm2.y()+=Foam::pow(stateRes[cellI].y(),2.0);
            vecResNorm2.z()+=Foam::pow(stateRes[cellI].z(),2.0);
        }
        totalResNorm2 += vecResNorm2.x() + vecResNorm2.y() + vecResNorm2.z();
    }
    
    forAll(adjReg_.volScalarStates,idxI)
    {
        const word stateName = adjReg_.volScalarStates[idxI];
        const word resName = stateName+"Res";                   
        const volScalarField& stateRes = db_.lookupObject<volScalarField>(resName); 
        
        scalar scalarResNorm2=0;
        forAll(stateRes,cellI)
        {
            scalarResNorm2+=Foam::pow(stateRes[cellI],2.0);
        }
        totalResNorm2 += scalarResNorm2;
    }
    
    forAll(adjRAS_.turbStates,idxI)
    {
        const word stateName = adjRAS_.turbStates[idxI];
        const word resName = stateName+"Res";  
        const volScalarField& stateRes = db_.lookupObject<volScalarField>(resName); 
        
        scalar scalarResNorm2=0;
        forAll(stateRes,cellI)
        {
            scalarResNorm2+=Foam::pow(stateRes[cellI],2.0);
        }
        totalResNorm2 += scalarResNorm2;
        turbResNorm2  += scalarResNorm2;
    }
    
    forAll(adjReg_.surfaceScalarStates,idxI)
    {
        const word stateName = adjReg_.surfaceScalarStates[idxI];
        const word resName = stateName+"Res";  
        const surfaceScalarField& stateRes = db_.lookupObject<surfaceScalarField>(resName); 
        
        forAll(stateRes,faceI)
        {
            phiResNorm2+=Foam::pow(stateRes[faceI],2.0);
    
        }
        forAll(stateRes.boundaryField(),patchI)
        {
            forAll(stateRes.boundaryField()[patchI],faceI)
            {
                scalar bPhiRes = stateRes.boundaryField()[patchI][faceI];
                phiResNorm2+=Foam::pow(bPhiRes,2.0);
            }
        }
        totalResNorm2 += phiResNorm2;
        
    }


    reduce(totalResNorm2,sumOp<scalar>());
    totalResNorm2=Foam::pow(totalResNorm2,0.5);

    reduce(turbResNorm2,sumOp<scalar>());
    turbResNorm2=Foam::pow(turbResNorm2,0.5);

    reduce(phiResNorm2,sumOp<scalar>());
    phiResNorm2=Foam::pow(phiResNorm2,0.5);

    if(mode=="total") return totalResNorm2;
    else if (mode=="turb") return turbResNorm2;
    else if (mode=="phi") return phiResNorm2;
    else FatalErrorIn("")<<"mode not valid"<< abort(FatalError);

    FatalErrorIn("")<<"mode not valid"<< abort(FatalError);
    return -10000.0;

}

void ReducedOrderModeling::getMatrixColNorms(Mat matIn,word matName)
{

    label nCols, nRows;
    MatGetSize(matIn,&nRows,&nCols);

    Vec colVec;
    VecCreate(PETSC_COMM_WORLD,&colVec);
    VecSetSizes(colVec,localSize_,PETSC_DECIDE);
    VecSetFromOptions(colVec);

    OFstream fOut(matName+"_ColInfo.txt");
    fOut<<"col                 L2Norm                Mean"<<endl;
    
    scalar norm2=0.0,mean=0.0;
    for(label n=0;n<nCols;n++)
    {
        VecZeroEntries(colVec);
        MatGetColumnVector(matIn,colVec,n);
        VecNorm(colVec,NORM_2,&norm2);
        VecSum(colVec,&mean);
        mean=mean*1.0/nRows;
        fOut<<n<<"  "<<norm2<<"  "<<mean<<endl;
    }

    return;
}

void ReducedOrderModeling::getPhiMatStateInfo(Mat matIn)
{
    label nRows,nCols;
    MatGetSize(matIn,&nRows,&nCols);

    label myID=Pstream::myProcNo();

    Vec vecOut;
    VecCreate(PETSC_COMM_WORLD,&vecOut);
    VecSetSizes(vecOut,nCols,PETSC_DECIDE);
    VecSetFromOptions(vecOut);

    forAll(adjReg_.volVectorStates,idxI)
    {
        const word stateName = adjReg_.volVectorStates[idxI];                
        const volVectorField& state = db_.lookupObject<volVectorField>(stateName); 
        
        scalar maxVolVectorMag=0.0;
        scalar maxVolVectorCellI=0.0;
        forAll(state,cellI)
        {
            if( mag(state[cellI])>maxVolVectorMag )
            {
                maxVolVectorMag=mag(state[cellI]);
                maxVolVectorCellI=cellI;
            }
        }

        VecZeroEntries(vecOut);

        for(label j=0;j<nCols;j++)
        {
            scalar val;
            vector vecMag=vector::zero;
            for (label i=0;i<3;i++)
            {
                label globalIdx = adjIdx_.getGlobalAdjointStateIndex(stateName,maxVolVectorCellI,i);
                MatGetValues(matIn,1,&globalIdx,1,&j,&val);
                vecMag[i]=val;
            }
            label globalRow=myID*nCols+j;
            scalar mag1=mag(vecMag);
            VecSetValue(vecOut,globalRow,mag1,INSERT_VALUES);
        }

        VecAssemblyBegin(vecOut);
        VecAssemblyEnd(vecOut);

        std::ostringstream nn("");
        nn<<stateName;
        std::string name1="SVDStateInfo_"+nn.str();
        adjIO_.writeVectorASCII(vecOut,name1);
    }

    
    return;
}

// end of namespace Foam
}

// ************************************************************************* //
