// meanToInstFields

#include "fvCFD.H"
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<<"Assigning Mean Fields to Instantaneous Fields"<<endl;

    argList::addOption
    (
        "varNames",
        "'(U p)'",
        "List of variable names to assign. Can be either a volScalar or volVector var."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // read options
    List<wordRe> varNames;
    if (args.optionFound("varNames"))
    {
        varNames = wordReList(args.optionLookup("varNames")());
    }
    else
    {
        Info<<"varNames not set! Exit."<<endl;
        return 0;
    }

    forAll(varNames,idxI)
    {
        word varName = varNames[idxI];
        Info<<"Assigning "<<varName<<endl;
        if (varName == "U")
        {
            volVectorField var
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            volVectorField varMean
            (
                IOobject
                (
                    "UMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            forAll(var,idxI)
            {
                var[idxI]=varMean[idxI];
            }
            var.correctBoundaryConditions();
            var.write();

        }
        else if (varName == "phi")
        {
            surfaceScalarField var
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            surfaceScalarField varMean
            (
                IOobject
                (
                    "phiMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            forAll(var,idxI)
            {
                var[idxI]=varMean[idxI];
            }
            forAll(var.boundaryField(),patchI)
            {
                forAll(var.boundaryField()[patchI],faceI)
                {
                    var.boundaryFieldRef()[patchI][faceI]=varMean.boundaryField()[patchI][faceI];
                }
            }
            var.write();
        }
        else if (varName=="nut")
        {

            volVectorField U
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            );
            #include "createPhi.H"
            singlePhaseTransportModel laminarTransport(U, phi);
            autoPtr<incompressible::turbulenceModel> turbulence
            (
                incompressible::turbulenceModel::New(U, phi, laminarTransport)
            );
            volScalarField var
            (
                IOobject
                (
                    varName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            volScalarField varMean
            (
                IOobject
                (
                    varName+"Mean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            forAll(var,idxI)
            {
                var[idxI]=varMean[idxI];
            }
            var.correctBoundaryConditions();
            var.write();
        }
        else
        {
            
            volScalarField var
            (
                IOobject
                (
                    varName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            volScalarField varMean
            (
                IOobject
                (
                    varName+"Mean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            forAll(var,idxI)
            {
                var[idxI]=varMean[idxI];
            }
            var.correctBoundaryConditions();
            var.write();
        }
        
    }
 
    return 0;
}


// ************************************************************************* //
