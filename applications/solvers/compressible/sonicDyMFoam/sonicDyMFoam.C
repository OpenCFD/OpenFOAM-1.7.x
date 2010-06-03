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
    sonicDyMFoam

Description
    Transient solver for trans-sonic/supersonic, laminar or turbulent flow
    of a compressible gas with mesh motion..

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    autoPtr<Foam::motionSolver> motionPtr = motionSolver::New(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "compressibleCourantNo.H"

        mesh.movePoints(motionPtr->newPoints());

        #include "rhoEqn.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U)
          + fvm::div(phi, U)
          + turbulence->divDevRhoReff(U)
        );

        solve(UEqn == -fvc::grad(p));

        #include "eEqn.H"


        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            U = UEqn.H()/UEqn.A();

            surfaceScalarField phid
            (
                "phid",
                fvc::interpolate(psi)
               *(
                    (fvc::interpolate(U) & mesh.Sf()) - fvc::meshPhi(rho, U)
                )
            );

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + fvm::div(phid, p)
                  - fvm::laplacian(rho/UEqn.A(), p)
                );

                pEqn.solve();

                phi = pEqn.flux();
            }

            #include "compressibleContinuityErrs.H"

            U -= fvc::grad(p)/UEqn.A();
            U.correctBoundaryConditions();
        }

        DpDt = fvc::DDt
        (
            surfaceScalarField
            (
                "phiU",
                phi/fvc::interpolate(rho) + fvc::meshPhi(rho, U)
            ),
            p
        );

        turbulence->correct();

        rho = psi*p;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
