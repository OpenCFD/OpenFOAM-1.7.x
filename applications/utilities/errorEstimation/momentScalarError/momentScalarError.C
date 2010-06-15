/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    momentScalarError

Description
    Estimates the error in the solution for a scalar transport equation in the
    standard form

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "linear.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    Info<< "\nEstimating error in scalar transport equation\n"
        << "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        IOobject THeader
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (THeader.headerOk() && Uheader.headerOk())
        {
            Info << "Reading T" << endl;
            volScalarField T(THeader, mesh);

            Info << "Reading U" << endl;
            volVectorField U(Uheader, mesh);

#           include "createPhi.H"

            volVectorField gradT = fvc::grad(T);

            volScalarField TE = 0.5*sqr(T);

            volScalarField L
            (
                IOobject
                (
                    "L",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("one", dimLength, 1.0)
            );

            L.internalField() =
                mesh.V()/fvc::surfaceSum(mesh.magSf())().internalField();

            // Divergence of the error in the T squared
            volScalarField momError
            (
                IOobject
                (
                    "momErrorL" + T.name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sqrt
                (
                    2.0*mag
                    (
                        (
                            fv::gaussConvectionScheme<scalar>
                            (
                                mesh,
                                phi,
                                tmp<surfaceInterpolationScheme<scalar> >
                                (
                                    new linear<scalar>(mesh)
                                )
                            ).fvcDiv(phi, TE)

                          - DT*
                            fv::gaussLaplacianScheme<scalar, scalar>(mesh)
                           .fvcLaplacian
                            (
                                TE
                            )
                          + DT*(gradT & gradT)
                        )*L/(mag(U) + DT/L)
                    )
                )
            );

            momError.boundaryField() = 0.0;
            momError.write();
        }
        else
        {
            Info<< "    No T or U" << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
