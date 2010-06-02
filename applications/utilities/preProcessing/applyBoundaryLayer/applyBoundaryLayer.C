/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    applyBoundaryLayer

Description
    Apply a simplified boundary-layer model to the velocity and
    turbulence fields based on the 1/7th power-law.

    The uniform boundary-layer thickness is either provided via the -ybl option
    or calculated as the average of the distance to the wall scaled with
    the thickness coefficient supplied via the option -Cbl.  If both options
    are provided -ybl is used.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("ybl", "scalar");
    argList::validOptions.insert("Cbl", "scalar");
    argList::validOptions.insert("writenut", "");

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading field U\n" << endl;
    volVectorField U
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

#   include "createPhi.H"

    Info<< "Calculating wall distance field" << endl;
    volScalarField y = wallDist(mesh).y();

    // Set the mean boundary-layer thickness
    dimensionedScalar ybl("ybl", dimLength, 0);

    if (args.optionFound("ybl"))
    {
        // If the boundary-layer thickness is provided use it
        ybl.value() = args.optionRead<scalar>("ybl");
    }
    else if (args.optionFound("Cbl"))
    {
        // Calculate boundary layer thickness as Cbl * mean distance to wall
        ybl.value() = gAverage(y) * args.optionRead<scalar>("Cbl");
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Neither option 'ybl' or 'Cbl' have been provided to calculate"
               " the boundary-layer thickness"
            << exit(FatalError);
    }

    Info<< "\nCreating boundary-layer for U of thickness "
        << ybl.value() << " m" << nl << endl;

    // Modify velocity by applying a 1/7th power law boundary-layer
    // u/U0 = (y/ybl)^(1/7)
    // assumes U0 is the same as the current cell velocity

    scalar yblv = ybl.value();
    forAll(U, celli)
    {
        if (y[celli] <= yblv)
        {
            U[celli] *= ::pow(y[celli]/yblv, (1.0/7.0));
        }
    }

    Info<< "Writing U" << endl;
    U.write();

    // Update/re-write phi
    phi = fvc::interpolate(U) & mesh.Sf();
    phi.write();

    // Set turbulence constants
    dimensionedScalar kappa("kappa", dimless, 0.41);
    dimensionedScalar Cmu("Cmu", dimless, 0.09);

    // Read and modify turbulence fields if present

    IOobject epsilonHeader
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject kHeader
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject nuTildaHeader
    (
        "nuTilda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // First calculate nut
    volScalarField nut
    (
        "nut",
        sqr(kappa*min(y, ybl))*::sqrt(2)*mag(dev(symm(fvc::grad(U))))
    );

    if (args.optionFound("writenut"))
    {
        Info<< "Writing nut" << endl;
        nut.write();
    }


    // Read and modify turbulence fields if present

    if (nuTildaHeader.headerOk())
    {
        Info<< "Reading field nuTilda\n" << endl;
        volScalarField nuTilda(nuTildaHeader, mesh);
        nuTilda = nut;
        nuTilda.correctBoundaryConditions();

        Info<< "Writing nuTilda\n" << endl;
        nuTilda.write();
    }

    if (kHeader.headerOk() && epsilonHeader.headerOk())
    {
        Info<< "Reading field k\n" << endl;
        volScalarField k(kHeader, mesh);

        Info<< "Reading field epsilon\n" << endl;
        volScalarField epsilon(epsilonHeader, mesh);

        scalar ck0 = ::pow(Cmu.value(), 0.25)*kappa.value();
        k = sqr(nut/(ck0*min(y, ybl)));
        k.correctBoundaryConditions();

        scalar ce0 = ::pow(Cmu.value(), 0.75)/kappa.value();
        epsilon = ce0*k*sqrt(k)/min(y, ybl);
        epsilon.correctBoundaryConditions();

        Info<< "Writing k\n" << endl;
        k.write();

        Info<< "Writing epsilon\n" << endl;
        epsilon.write();
    }

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
