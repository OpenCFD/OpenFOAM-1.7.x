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
    rhoSonicFoam

Description
    Density-based compressible flow solver.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readThermodynamicProperties.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        surfaceScalarField phiv
        (
            IOobject
            (
                "phiv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rhoU)/fvc::interpolate(rho) & mesh.Sf()
        );

        scalar CoNum = max
        (
            mesh.surfaceInterpolation::deltaCoeffs()
           *mag(phiv)/mesh.magSf()
        ).value()*runTime.deltaT().value();

        Info<< "\nMax Courant Number = " << CoNum << endl;

        solve
        (
            fvm::ddt(rho)
          + fvm::div(phiv, rho)
        );

        p = rho/psi;

        solve
        (
            fvm::ddt(rhoU)
          + fvm::div(phiv, rhoU)
         == 
          - fvc::grad(p)
        );

        U == rhoU/rho;

        surfaceScalarField phiv2
        (
            IOobject
            (
                "phiv2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rhoU)/fvc::interpolate(rho) & mesh.Sf()
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvm::div(phiv, rhoE)
         ==
          - fvc::div(phiv2, p)
        );

        T = (rhoE - 0.5*rho*magSqr(rhoU/rho))/Cv/rho;
        psi = 1.0/(R*T);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
