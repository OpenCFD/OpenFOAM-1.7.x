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
    dieselEngineFoam

Description
    Solver for diesel engine spray and combustion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "engineTime.H"
#include "engineMesh.H"
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "spray.H"
#include "psiChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"
#include "Switch.H"
#include "OFstream.H"
#include "volPointInterpolation.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createEngineTime.H"
    #include "createEngineMesh.H"
    #include "createFields.H"
    #include "readGravitationalAcceleration.H"
    #include "readCombustionProperties.H"
    #include "createSpray.H"
    #include "initContinuityErrs.H"
    #include "readEngineTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "startSummary.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readPISOControls.H"
        #include "readEngineTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Crank angle = " << runTime.theta() << " CA-deg" << endl;

        mesh.move();
        const_cast<volPointInterpolation&>
        (
            volPointInterpolation::New(mesh)
        ).updateMesh();

        dieselSpray.evolve();

        Info<< "Solving chemistry" << endl;

        chemistry.solve
        (
            runTime.value() - runTime.deltaT().value(),
            runTime.deltaT().value()
        );

        // turbulent time scale
        {
            volScalarField tk =
                Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
            volScalarField tc = chemistry.tc();

            //Chalmers PaSR model
            kappa = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
        }

        chemistrySh = kappa*chemistry.Sh()();

        #include "rhoEqn.H"
        #include "UEqn.H"

        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
            #include "YEqn.H"
            #include "hsEqn.H"

            // --- PISO loop
            for (int corr=1; corr<=nCorr; corr++)
            {
                #include "pEqn.H"
            }
        }

        turbulence->correct();

        #include "logSummary.H"
        #include "spraySummary.H"

        rho = thermo.rho();

        if (runTime.write())
        {
            chemistry.dQ()().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
