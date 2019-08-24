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
    AdjointDerivative& adjDev
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
    mfStep                    = readOptionOrDefault<scalar>(romDict_,"mfStep",1e-8);

    // print all the parameters to screen    
    Info<<"ROM Parameters"<<romParameters_<<endl;

    localSize_=adjIdx_.nLocalAdjointStates;
    nFFDs_=adjIO_.nFFDPoints;

}

ReducedOrderModeling::~ReducedOrderModeling()
{
}

void ReducedOrderModeling::initializeOffline()
{
    // save the unperturb residual statistics as reference
    // we will verify against this reference to ensure consistent residuals after perturbing
    // and resetting states
    adjDev_.calcFlowResidualStatistics("set");

    this->initializeSnapshotMat();

}

void ReducedOrderModeling::initializeOnline()
{
    
    this->initializedRdWMatReduced();
    this->initializedRdFFDMatReduced();
    this->initializeSVDPhiMat();

    label nProcs = Pstream::nProcs();
    std::ostringstream np("");
    np<<nProcs;
    std::string fNamedRdWReduced="dRdWReduced_"+np.str();
    std::string fNamedRdFFDReduced="dRdFFDReduced_"+np.str();
    std::string fNamePhi="svdPhiMat_"+np.str();

    adjIO_.readMatrixBinary(dRdWReduced_,fNamedRdWReduced);
    adjIO_.readMatrixBinary(dRdFFDReduced_,fNamedRdFFDReduced);
    adjIO_.readMatrixBinary(svdPhiMat_,fNamePhi);
    
    VecCreate(PETSC_COMM_WORLD,&deltaFFDVec_);
    VecSetSizes(deltaFFDVec_,PETSC_DECIDE,nFFDs_);
    VecSetFromOptions(deltaFFDVec_);

    label Istart,Iend;
    VecGetOwnershipRange(deltaFFDVec_,&Istart,&Iend);

    for(label i=Istart;i<Iend;i++)
    {
        scalar val=deltaFFD[i];
        VecSetValue(deltaFFDVec_,i,val,INSERT_VALUES);
    }
    VecAssemblyBegin(deltaFFDVec_);
    VecAssemblyEnd(deltaFFDVec_);
    

}

void ReducedOrderModeling::initializeSnapshotMat()
{
    Info<<"Initializing the Snapshot matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // now initialize the memory for the jacobian itself
    // create snapshotMat_
    MatCreate(PETSC_COMM_WORLD,&snapshotMat_);
    MatSetSizes(snapshotMat_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(snapshotMat_);
    MatMPIAIJSetPreallocation(snapshotMat_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(snapshotMat_,nSamples,NULL);
    MatSetUp(snapshotMat_);
    Info<<"Snapshot matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    return;
}

void ReducedOrderModeling::initializedRdFFDMat()
{
    //******************************** Compute dRdFFD *****************************//
    Info<<"Initializing the dRdFFD matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // create dRdFFD_
    MatCreate(PETSC_COMM_WORLD,&dRdFFD_);
    MatSetSizes(dRdFFD_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,adjIO_.nFFDPoints);
    MatSetFromOptions(dRdFFD_);
    MatMPIAIJSetPreallocation(dRdFFD_,adjIO_.nFFDPoints,NULL,adjIO_.nFFDPoints,NULL);
    MatSeqAIJSetPreallocation(dRdFFD_,adjIO_.nFFDPoints,NULL);
    MatSetOption(dRdFFD_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(dRdFFD_);
    Info<<"dRdFFD matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    return;
}

void ReducedOrderModeling::initializeSVDPhiMat()
{
    Info<<"Initializing the Phi matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // create svdPhiMat_
    MatCreate(PETSC_COMM_WORLD,&svdPhiMat_);
    MatSetSizes(svdPhiMat_,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(svdPhiMat_);
    MatMPIAIJSetPreallocation(svdPhiMat_,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(svdPhiMat_,nSamples,NULL);
    MatSetOption(svdPhiMat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(svdPhiMat_);
    Info<<"Phi matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;

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
    Info<<"Setting the snapshot matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;

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
                    MatSetValues(snapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
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
                MatSetValues(snapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
            }          
        }
    
        // perturb turbStates
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
                MatSetValues(snapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
            }          
    
        }
        // perturb surfaceScalarStates
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
                MatSetValues(snapshotMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
            }
        }

    }

    MatAssemblyBegin(snapshotMat_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(snapshotMat_,MAT_FINAL_ASSEMBLY);

    //adjIO_.writeMatrixASCII(snapshotMat_,"snapshotMat");

    Info<<"The snapshot matrix is set. "<<runTime_.elapsedClockTime()<<" s"<<endl;
}

void ReducedOrderModeling::solveOffline()
{
    PetscInt       nsv,maxit,nconv,its;
    PetscScalar    tol,sigma,val;
    SVDType        type;

    adjDev_.calcFlowResidualStatistics("print");

    this->setSnapshotMat();

    MatCreateVecs(snapshotMat_,&vVec_,&uVec_);

    this->initializeSVDPhiMat();

    PetscInt Istart,Iend;
    MatGetOwnershipRange(svdPhiMat_,&Istart,&Iend);


    Info<<"Solving the SVD..."<<endl;
    SVD svd;
    SVDCreate(PETSC_COMM_WORLD,&svd);
    SVDSetOperator(svd,snapshotMat_);
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
        SVDGetSingularTriplet(svd,i,&sigma,uVec_,vVec_);
        Info<<"  "<<sigma<<endl;
        PetscScalar *uVecArray;
        VecGetArray(uVec_,&uVecArray);
        for(label j=Istart;j<Iend;j++)
        {
            label relIdx = j-Istart;
            val = uVecArray[relIdx];
            MatSetValues(svdPhiMat_,1,&j,1,&i,&val,INSERT_VALUES);
        }
        VecRestoreArray(uVec_,&uVecArray);
    }
    MatAssemblyBegin(svdPhiMat_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(svdPhiMat_,MAT_FINAL_ASSEMBLY);

    Info<<"Writing the phiMat"<<endl;
    label nProcs = Pstream::nProcs();
    std::ostringstream np("");
    np<<nProcs;
    std::string fNamePhi="svdPhiMat_"+np.str();

    adjIO_.writeMatrixBinary(svdPhiMat_,fNamePhi);
    //adjIO_.writeMatrixASCII(svdPhiMat_,"svdPhiMat");

    MatDestroy(&snapshotMat_);
    SVDDestroy(&svd);

    Info<<"Calculating dRdW and dRdFFD at time = "<<runTime_.timeName()<<endl;

    
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
        Mat dRdWsvdPhiMat;    
        Info<< "Computing phiT*dRdW*phi" << endl;
        // now compute the matmat mult
        MatMatMult(dRdW_,svdPhiMat_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdWsvdPhiMat);
        MatTransposeMatMult(svdPhiMat_,dRdWsvdPhiMat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdWReduced_);
    
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
        MatTransposeMatMult(svdPhiMat_,dRdFFD_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdFFDReduced_);
    }

    Info<<"Writing the reduced matrices...."<<endl;

    std::string fNamedRdWReduced="dRdWReduced_"+np.str();
    adjIO_.writeMatrixBinary(dRdWReduced_,fNamedRdWReduced);
    adjIO_.writeMatrixASCII(dRdWReduced_,fNamedRdWReduced);

    std::string fNamedRdFFDReduced="dRdFFDReduced_"+np.str();
    adjIO_.writeMatrixBinary(dRdFFDReduced_,fNamedRdFFDReduced);
    adjIO_.writeMatrixASCII(dRdFFDReduced_,fNamedRdFFDReduced);

    Info<<"Writing the reduced matrices.... Done!"<<endl;

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
    MatGetColumnVector(svdPhiMat_,colVec,n);

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

void ReducedOrderModeling::calcReducedMatsMF()
{

    // ********************** Ar *********************
    // first compute dRdW*phi
    Info<<"Computing dRdW*Phi..."<<endl;
    Mat dRdWPhi;
    MatCreate(PETSC_COMM_WORLD,&dRdWPhi);
    MatSetSizes(dRdWPhi,localSize_,PETSC_DECIDE,PETSC_DETERMINE,nSamples);
    MatSetFromOptions(dRdWPhi);
    MatMPIAIJSetPreallocation(dRdWPhi,nSamples,NULL,nSamples,NULL);
    MatSeqAIJSetPreallocation(dRdWPhi,nSamples,NULL);
    MatSetUp(dRdWPhi);

    // do matrix-free for dRdW*Phi = [ R(w+v*mfSetp) - R(w) ] / mfStep
    label isRef=1, isPC=0;
    adjDev_.copyStates("Ref2Var");
    adjDev_.calcResiduals(isRef,isPC);
    adjRAS_.calcTurbResiduals(isRef,isPC);
    for(label nn=0;nn<nSamples;nn++)
    {
        this->perturbStatesMF(nn);

        adjDev_.calcResPartDeriv(mfStep,isPC);

        // reset perturbation
        adjDev_.copyStates("Ref2Var");

        this->setdRdWPhiMat(dRdWPhi,nn);

    }

    MatAssemblyBegin(dRdWPhi,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(dRdWPhi,MAT_FINAL_ASSEMBLY);

    adjDev_.calcFlowResidualStatistics("verify");

    // now compute phiT*dRdW*Phi
    //Mat phiT;
    //MatCreateTranspose(svdPhiMat_,&phiT);
    void initializedRdWMatReduced();
    Info<< "Computing phiT*dRdW*phi" << endl;
    MatTransposeMatMult(svdPhiMat_,dRdWPhi,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdWReduced_);

    // ********************** Br *********************
    Info<<"Calculating dRdFFD... "<<endl;
    this->initializedRdFFDMat();
    adjDev_.calcdRdFFD(dRdFFD_);
    Info<< "Computing phiT*dRdFFD" << endl;
    void initializedRdFFDMatReduced();
    MatTransposeMatMult(svdPhiMat_,dRdFFD_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dRdFFDReduced_);

    MatDestroy(&dRdFFD_);
    
    return;

}

void ReducedOrderModeling::setNewField(Vec deltaWVec)
{
    const objectRegistry& db_(mesh_.thisDb());

    label Istart,Iend;
    label cellI,comp,faceI;
    VecGetOwnershipRange(deltaWVec,&Istart,&Iend);

    PetscScalar *deltaWArray;
    VecGetArray(deltaWVec,&deltaWArray);

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
            state[cellI][comp] += deltaWArray[relIdx];
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
            state[cellI] += deltaWArray[relIdx];
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
                state[faceI] += deltaWArray[relIdx];
            }
            else
            {
                label relIdx=faceI-adjIdx_.nLocalInternalFaces;
                label patchIdx=adjIdx_.bFacePatchI[relIdx];
                label faceIdx=adjIdx_.bFaceFaceI[relIdx];
                state.boundaryFieldRef()[patchIdx][faceIdx] += deltaWArray[relIdx];
            }
        }
        else
        {
            FatalErrorIn("")<<"stateType not known"<< abort(FatalError);
        }
    }
    
    VecRestoreArray(deltaWVec,&deltaWArray);

    // need to update BCs and intermediate variable
    adjDev_.updateStateVariableBCs();
    adjDev_.updateIntermediateVariables();
}

void ReducedOrderModeling::writeNewField()
{
    const objectRegistry& db_(mesh_.thisDb());
    // write
    forAll(adjReg_.volVectorStates,idxI)                                           
    {        
        // create state and stateRef
        makeState(volVectorStates,volVectorField,adjReg_); 
        word oldName=state.name();
        word newName=state.name()+"ROM";
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
        word newName=state.name()+"ROM";
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
        word newName=state.name()+"ROM";
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
    nut.rename("nutROM");
    nut.write();
    nut.rename("nut");

    forAll(adjReg_.surfaceScalarStates,idxI)
    {
        // create state and stateRef
        makeState(surfaceScalarStates,surfaceScalarField,adjReg_);
        word oldName=state.name();
        word newName=state.name()+"ROM";
        state.rename(newName);
        state.write();       
        state.rename(oldName); 
    }

}

void ReducedOrderModeling::solveOnline()
{
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
 
    VecCreate(PETSC_COMM_WORLD,&deltaWVec);
    VecSetSizes(deltaWVec,localSize_,PETSC_DECIDE);
    VecSetFromOptions(deltaWVec);

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

    MatMult(svdPhiMat_,deltaWVecReduced,deltaWVec);

    //adjIO_.writeVectorASCII(deltaWVecReduced,"deltaWVecReduced");
    //adjIO_.writeVectorASCII(deltaWVec,"deltaWVec");

    this->setNewField(deltaWVec);
    adjObj_.writeObjFuncValues();
    this->writeNewField();

}

// end of namespace Foam
}

// ************************************************************************* //
