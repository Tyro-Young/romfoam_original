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

    Mat snapshotMat_,dRdW_,dRdFFD_,svdPhiMat_, dRdWReduced_, dRdFFDReduced_;
    Vec uVec_,vVec_,deltaFFDVec_;

    label nFFDs_;
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

    label nSamples;
    scalarList deltaFFD;
    word svdType;
    scalar svdTol;
    label svdMaxIts;
    label svdRequestedN;
    label useMF;
    scalar mfStep;

    void initializeSnapshotMat();
    void initializedRdFFDMat();
    void initializeSVDPhiMat();
    void initializedRdWMatReduced();
    void initializedRdFFDMatReduced();
    
    void initializeOffline();
    void solveOffline();
    void setSnapshotMat();
    void calcReducedMatsMF();
    void perturbStatesMF(label n);
    void setdRdWPhiMat(Mat matIn,label n);


    void initializeOnline();
    void solveOnline();

    void setNewField(Vec deltaWVec);
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
