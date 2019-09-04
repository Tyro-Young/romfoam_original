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

#ifndef ReducedOrderModeling_H
#define ReducedOrderModeling_H

#include "fvOptions.H"
#include "surfaceFields.H"
#include <slepcsvd.h>
#include "AdjointIO.H"
#include "AdjointSolverRegistry.H"
#include "AdjointRASModel.H"
#include "AdjointIndexing.H"
#include "AdjointJacobianConnectivity.H"
#include "AdjointObjectiveFunction.H"
#include "AdjointDerivative.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ReducedOrderModeling Declaration
\*---------------------------------------------------------------------------*/
class ReducedOrderModeling
{

private:

    // Disallow default bitwise copy construct and assignment
    ReducedOrderModeling(const ReducedOrderModeling&);
    void operator=(const ReducedOrderModeling&);

    Foam::Time& runTime_;
    fvMesh& mesh_;

    AdjointIO& adjIO_;
    AdjointSolverRegistry& adjReg_;
    AdjointRASModel& adjRAS_;
    AdjointIndexing& adjIdx_;
    AdjointJacobianConnectivity& adjCon_;
    AdjointObjectiveFunction& adjObj_;
    AdjointDerivative& adjDev_;

protected:

    /// ROM dict
    IOdictionary romDict_;

    /// A temp dict to store all the ROM parameters and print them to screen later 
    IOdictionary romParameters_; 

    /// snapshot matrix for all state variables
    Mat snapshotMat_;
    
    /// partial R / partial W matrix
    Mat dRdW_;
    
    /// partial R / partial FFD matrix
    Mat dRdFFD_;
    
    /// the left matrix resulted from the SVD computation
    Mat svdPhiMat_;
    
    /// reduced dRdW matrix
    Mat dRdWReduced_;
    
    /// reduced dRdFFD matrix
    Mat dRdFFDReduced_;
    
    /// left vector for SVD uVec*PhiMat*vVec with size same as the row size for PhiMat
    Vec uVec_;
    
    /// right vector for SVD, uVec*PhiMat*vVec with size same as the column size for PhiMat
    Vec vVec_;
    
    /// vector to store the delta FFD displacement for prediction
    Vec deltaFFDVec_;

    /// number of FFD points or design variables
    label nFFDs_;
    
    /// size of local state variables
    label localSize_;

public:

    // Constructors

    //- Construct from components
    ReducedOrderModeling
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
    );


    //- Destructor
    virtual ~ReducedOrderModeling();

    // we make these variable public
    /// number of snapshot samples
    label nSamples;

    /// list of delta FFD displacement for prediction
    scalarList deltaFFD;

    /// type of SVD computation
    word svdType;

    /// SVD convergence tolerance
    scalar svdTol;

    /// SVD max iteration
    label svdMaxIts;

    /// required SVD component
    label svdRequestedN;

    /// whether to use matrix free method for ROM
    label useMF;

    /// finite-difference step size for the matrix free method
    scalar mfStep;

    /// whether to use the debug mode where we output some vecs and matc
    label debugMode;

    /// initialize snapshot matrix
    void initializeSnapshotMat();

    /// initialize dRdFFD
    void initializedRdFFDMat();

    /// initialize svdPhiMat
    void initializeSVDPhiMat();

    /// initialize dRdWReduced
    void initializedRdWMatReduced();

    /// initialize dRdFFDReduced
    void initializedRdFFDMatReduced();
    
    /// initialize offlineROM
    void initializeOffline();

    /// solver offline ROM
    void solveOffline();

    /// assign values for snapshotMat
    void setSnapshotMat();

    /// compute reduced matrices using matrix-free method
    void calcReducedMatsMF();

    /// perturb state variables using matrix-free method
    void perturbStatesMF(label n);

    /// set values for the matrix-matrix product dRdW*Phi 
    void setdRdWPhiMat(Mat matIn,label n);

    /// initialize onlineROM
    void initializeOnline();

    /// solve onlineROM
    void solveOnline();

    /// set the new fields (U p nuTilda etc) for ROM
    void setNewField(Vec deltaWVec);

    /// write the new fields as UROM pROM etc.
    void writeNewField();

    /// read an option from dict and add them to romParameters_
    template<class type>
    type readOptionOrDefault
    (
        const dictionary& dict,
        word option,
        type defaultValue
    );



};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class type>
type ReducedOrderModeling::readOptionOrDefault
(
    const dictionary& dict,
    word option,
    type defaultValue
)
{
    
    type value = dict.lookupOrDefault<type>(option,defaultValue); 
    if (value == defaultValue)
    {
        // notify the user the default value is used.
        word optionPrint = option+"(D)"; 
        romParameters_.add(optionPrint,value);
    }
    else
    {
        romParameters_.add(option,value);
    }
    return value;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
