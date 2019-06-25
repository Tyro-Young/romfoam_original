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
    timeSamples               = readOptionOrDefault<labelList>(romDict_,"timeSamples",{});

    // print all the parameters to screen    
    Info<<"ROM Parameters"<<romParameters_<<endl;

}

ReducedOrderModeling::~ReducedOrderModeling()
{
}

void ReducedOrderModeling::initialize()
{
    this->initializeMats();
    this->initializeVecs();

}

void ReducedOrderModeling::initializeMats()
{
    Info<<"Initializing the Snapshot matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // now initialize the memory for the jacobian itself
    label localSize = adjIdx_.nLocalAdjointStates;
    label nInstances = timeSamples.size();
    // create AMat_
    MatCreate(PETSC_COMM_WORLD,&AMat_);
    MatSetSizes(AMat_,localSize,PETSC_DECIDE,PETSC_DETERMINE,nInstances);
    MatSetFromOptions(AMat_);
    MatMPIAIJSetPreallocation(AMat_,nInstances,NULL,nInstances,NULL);
    MatSeqAIJSetPreallocation(AMat_,nInstances,NULL);
    MatSetOption(AMat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(AMat_);
    Info<<"Snapshot matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    Info<<"Initializing the Phi matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // create phiMat_
    MatCreate(PETSC_COMM_WORLD,&phiMat_);
    MatSetSizes(phiMat_,localSize,PETSC_DECIDE,PETSC_DETERMINE,nInstances);
    MatSetFromOptions(phiMat_);
    MatMPIAIJSetPreallocation(phiMat_,nInstances,NULL,nInstances,NULL);
    MatSeqAIJSetPreallocation(phiMat_,nInstances,NULL);
    MatSetOption(phiMat_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(phiMat_);
    Info<<"Phi matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    return;
}

void ReducedOrderModeling::initializeVecs()
{
    Info<<"Initializing vectors. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    MatCreateVecs(AMat_,&vVec_,&uVec_);

    label nInstances = timeSamples.size();
    VecCreate(PETSC_COMM_WORLD,&sigmaVec_);
    VecSetSizes(sigmaVec_,PETSC_DECIDE,nInstances);
    VecSetFromOptions(sigmaVec_);

    Info<<"Vectors created. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    return;
}

void ReducedOrderModeling::setSnapshotMat()
{
    Info<<"Setting the snapshot matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    forAll(timeSamples,idxI)
    {
        Info<< "Adding variables for time = " << name(timeSamples[idxI]) << endl;

        volVectorField U
        (
            IOobject
            (
                "U",
                name(timeSamples[idxI]),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );
        forAll(U,cellI)
        {
            for(label i=0;i<3;i++)
            {
                PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex("U",cellI,i);
                PetscInt colI = idxI;
                PetscScalar val = U[cellI][i];
                MatSetValues(AMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
            }
        }

        volScalarField p
        (
            IOobject
            (
                "p",
                name(timeSamples[idxI]),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        forAll(p,cellI)
        {
            PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex("p",cellI);
            PetscInt colI = idxI;
            PetscScalar val = p[cellI];
            MatSetValues(AMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
        }

        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                name(timeSamples[idxI]),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        forAll(mesh_.faces(), faceI)
        {
            PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex("phi",faceI);
            PetscInt colI = idxI;
            PetscScalar val;
            if (faceI < adjIdx_.nLocalInternalFaces)
            {
                val = phi[faceI];
            }
            else
            {
                label relIdx=faceI-adjIdx_.nLocalInternalFaces;
                label patchIdx=adjIdx_.bFacePatchI[relIdx];
                label faceIdx=adjIdx_.bFaceFaceI[relIdx];
                val = phi.boundaryField()[patchIdx][faceIdx];
            } 
            MatSetValues(AMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
  
        }


        volScalarField nuTilda
        (
            IOobject
            (
                "nuTilda",
                name(timeSamples[idxI]),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        forAll(nuTilda,cellI)
        {
            PetscInt rowI = adjIdx_.getGlobalAdjointStateIndex("nuTilda",cellI);
            PetscInt colI = idxI;
            PetscScalar val = nuTilda[cellI];
            MatSetValues(AMat_,1,&rowI,1,&colI,&val,INSERT_VALUES);
        }

    }

    Info<< "Begin assemblying Amat" << endl;
    MatAssemblyBegin(AMat_,MAT_FINAL_ASSEMBLY);
    Info<< "End assemblying Amat" << endl;
    MatAssemblyEnd(AMat_,MAT_FINAL_ASSEMBLY);

    Info<<"The snapshot matrix is set. "<<runTime_.elapsedClockTime()<<" s"<<endl;
}

void ReducedOrderModeling::solve()
{
    PetscInt       nsv,maxit,nconv,its;
    PetscScalar    tol,sigma,val;
    SVDType        type;

    this->setSnapshotMat();

    PetscInt Istart,Iend;
    MatGetOwnershipRange(phiMat_,&Istart,&Iend);

    SVDCreate(PETSC_COMM_WORLD,&svd_);
    SVDSetOperator(svd_,AMat_);
    SVDSetFromOptions(svd_);

    SVDSolve(svd_);

    // Output SVD solution diagnostics
    SVDGetIterationNumber(svd_,&its);
    PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);

    SVDGetType(svd_,&type);
    PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);

    SVDGetDimensions(svd_,&nsv,NULL,NULL);
    PetscPrintf(PETSC_COMM_WORLD," Number of requested singular values: %D\n",nsv);

    SVDGetTolerances(svd_,&tol,&maxit);
    PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);

    SVDGetConverged(svd_,&nconv);
    PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate singular triplets: %D\n\n",nconv);

    if (nconv>0) 
    {
        // Display singular values and relative errors
        PetscPrintf(PETSC_COMM_WORLD,
         "          sigma           relative error\n"
         "  --------------------- ------------------\n");
        for (label i=0;i<nconv;i++) 
        {
            // Get converged singular triplets: i-th singular value is stored in sigma
            SVDGetSingularTriplet(svd_,i,&sigma,uVec_,vVec_);

            Info<<sigma<<endl;

            VecSetValues(sigmaVec_,1,&i,&sigma,INSERT_VALUES);

            PetscScalar *uVecArray;
            VecGetArray(uVec_,&uVecArray);

            for(label j=Istart;j<Iend;j++)
            {
                label relIdx = j-Istart;
                val = uVecArray[relIdx];
                MatSetValues(phiMat_,1,&j,1,&i,&val,INSERT_VALUES);
            }

            VecRestoreArray(uVec_,&uVecArray);
        }

        VecAssemblyBegin(sigmaVec_);
        VecAssemblyEnd(sigmaVec_);

        MatAssemblyBegin(phiMat_,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(phiMat_,MAT_FINAL_ASSEMBLY);
    }

    //******************************** Compute dRdW *****************************//
    // compute preallocation vecs for dRdW 
    adjCon_.setupdRdWCon(1);
    adjCon_.initializedRdWCon();
    adjCon_.setupdRdWCon(0);
    // Read in the precomputed coloring
    adjCon_.readdRdWColoring();
    
    label transposed=0,isPC=0;
    adjDev_.initializedRdW(&dRdW_,transposed);

    adjDev_.calcdRdW(dRdW_,transposed,isPC);

    adjCon_.deletedRdWCon();

    //******************************** Compute Ar *****************************//
    Mat AMatR, phiMatT, phiMatTAMat,BMatR;

    Info<< "Transposing Phi to PhiT" << endl;
    MatTranspose(phiMat_,MAT_INITIAL_MATRIX,&phiMatT);
    MatAssemblyBegin(phiMatT,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(phiMatT,MAT_FINAL_ASSEMBLY);

    // Compute first part of Ar (PhiT*A)
    Info<< "Computing PhiT*A" << endl;
    MatMatMult(phiMatT,dRdW_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&phiMatTAMat);

    // Compute second part of Ar (PhiT*A*Phi)
    Info<< "Computing PhiTA*Phi" << endl;
    MatMatMult(phiMatTAMat,phiMat_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&AMatR);

    adjIO_.writeMatrixBinary(AMatR,"AMatR");
    adjIO_.writeMatrixASCII(AMatR,"AMatR");

    MatDestroy(&dRdW_);

    //******************************** Compute dRdFFD *****************************//
    Info<<"Initializing the dRdFFD matrix. "<<runTime_.elapsedClockTime()<<" s"<<endl;
    // create dRdFFD_
    label localSize = adjIdx_.nLocalAdjointStates;
    MatCreate(PETSC_COMM_WORLD,&dRdFFD_);
    MatSetSizes(dRdFFD_,localSize,PETSC_DECIDE,PETSC_DETERMINE,adjIO_.nFFDPoints);
    MatSetFromOptions(dRdFFD_);
    MatMPIAIJSetPreallocation(dRdFFD_,adjIO_.nFFDPoints,NULL,adjIO_.nFFDPoints,NULL);
    MatSeqAIJSetPreallocation(dRdFFD_,adjIO_.nFFDPoints,NULL);
    MatSetOption(dRdFFD_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(dRdFFD_);
    Info<<"dRdFFD matrix Created. "<<runTime_.elapsedClockTime()<<" s"<<endl;

    adjDev_.calcdRdFFD(dRdFFD_);

    //******************************** Compute Br *****************************//
    Info<< "Computing PhiT*B" << endl;
    MatMatMult(phiMatT,dRdFFD_,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&BMatR);

    adjIO_.writeMatrixBinary(BMatR,"BMatR");
    adjIO_.writeMatrixASCII(BMatR,"BMatR");


}


// end of namespace Foam
}

// ************************************************************************* //
