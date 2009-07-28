/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    stressComponents

Description
    Calculates and writes the scalar fields of the six components of the stress
    tensor sigma for each time.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check U exists
        if (Uheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

#           include "createPhi.H"

            singlePhaseTransportModel laminarTransport(U, phi);

            volSymmTensorField sigma
            (
                IOobject
                (
                    "sigma",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                laminarTransport.nu()*2*dev(symm(fvc::grad(U)))
            );


            volScalarField sigmaxx
            (
                IOobject
                (
                    "sigmaxx",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                sigma.component(tensor::XX)
            );
            sigmaxx.write();

            volScalarField sigmayy
            (
                IOobject
                (
                    "sigmayy",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                sigma.component(tensor::YY)
            );
            sigmayy.write();

            volScalarField sigmazz
            (
                IOobject
                (
                    "sigmazz",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                sigma.component(tensor::ZZ)
            );
            sigmazz.write();

            volScalarField sigmaxy
            (
                IOobject
                (
                    "sigmaxy",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                sigma.component(tensor::XY)
            );
            sigmaxy.write();

            volScalarField sigmaxz
            (
                IOobject
                (
                    "sigmaxz",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                sigma.component(tensor::XZ)
            );
            sigmaxz.write();

            volScalarField sigmayz
            (
                IOobject
                (
                    "sigmayz",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                sigma.component(tensor::YZ)
            );
            sigmayz.write();

            volVectorField Ub
            (
                IOobject
                (
                    "Ub",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                U,
                zeroGradientFvPatchVectorField::typeName
            );
            Ub.correctBoundaryConditions();
            Ub.write();

            volScalarField sigmaUn
            (
                IOobject
                (
                    "sigmaUn",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                0.0*sigma.component(tensor::YZ)
            );

            forAll(sigmaUn.boundaryField(), patchI)
            {
                sigmaUn.boundaryField()[patchI] =
                (
                    mesh.boundary()[patchI].nf()
                  & sigma.boundaryField()[patchI]
                )().component(vector::X);
            }

            sigmaUn.write();
        }
        else
        {
            Info<< "    No U" << endl;
        }

        Info<< endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
