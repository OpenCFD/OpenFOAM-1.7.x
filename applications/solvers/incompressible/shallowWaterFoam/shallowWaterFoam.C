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
    shallowWaterFoam

Description
    Transient solver for inviscid shallow-water equations with rotation.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"

        for (int ucorr=0; ucorr<nOuterCorr; ucorr++)
        {
            surfaceScalarField phiv("phiv", phi/fvc::interpolate(h));

            fvVectorMatrix hUEqn
            (
                fvm::ddt(hU)
              + fvm::div(phiv, hU)
            );

            hUEqn.relax();

            if (momentumPredictor)
            {
                if (rotating)
                {
                    solve(hUEqn + (F ^ hU) == -magg*h*fvc::grad(h + h0));
                }
                else
                {
                    solve(hUEqn == -magg*h*fvc::grad(h + h0));
                }

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                    hU.correctBoundaryConditions();
                }
            }

            // --- PISO loop
            for (int corr=0; corr<nCorr; corr++)
            {
                surfaceScalarField hf = fvc::interpolate(h);
                volScalarField rUA = 1.0/hUEqn.A();
                surfaceScalarField ghrUAf = magg*fvc::interpolate(h*rUA);

                surfaceScalarField phih0 = ghrUAf*mesh.magSf()*fvc::snGrad(h0);

                if (rotating)
                {
                    hU = rUA*(hUEqn.H() - (F ^ hU));
                }
                else
                {
                    hU = rUA*hUEqn.H();
                }

                phi = (fvc::interpolate(hU) & mesh.Sf())
                    + fvc::ddtPhiCorr(rUA, h, hU, phi)
                    - phih0;

                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    fvScalarMatrix hEqn
                    (
                        fvm::ddt(h)
                      + fvc::div(phi)
                      - fvm::laplacian(ghrUAf, h)
                    );

                    if (ucorr < nOuterCorr-1 || corr < nCorr-1)
                    {
                        hEqn.solve();
                    }
                    else
                    {
                        hEqn.solve(mesh.solver(h.name() + "Final"));
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi += hEqn.flux();
                    }
                }

                hU -= rUA*h*magg*fvc::grad(h + h0);

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                }

                hU.correctBoundaryConditions();
            }
        }

        U == hU/h;
        hTotal == h + h0;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
