/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    porousExplicitSourceReactingParcelFoam

Description
    Transient PISO solver for compressible, laminar or turbulent flow with
    reacting multiphase Lagrangian parcels for porous media, including explicit
    sources for mass, momentum and energy

    The solver includes:
    - reacting multiphase parcel cloud
    - porous media
    - mass, momentum and energy sources
    - polynomial based, incompressible thermodynamics (f(T))

    Note: ddtPhiCorr not used here when porous zones are active
    - not well defined for porous calculations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hReactionThermo.H"
#include "turbulenceModel.H"
#include "BasicReactingMultiphaseCloud.H"
#include "rhoChemistryModel.H"
#include "chemistrySolver.H"
#include "radiationModel.H"
#include "porousZones.H"
#include "timeActivatedExplicitSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readChemistryProperties.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createRadiationModel.H"
    #include "createClouds.H"
    #include "createExplicitSources.H"
    #include "createPorousZones.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readPISOControls.H"
        #include "readAdditionalSolutionControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        #include "chemistry.H"
        #include "rhoEqn.H"
        #include "UEqn.H"
        #include "YEqn.H"
        #include "hsEqn.H"

        // --- PISO loop
        for (int corr=0; corr<nCorr; corr++)
        {
            #include "pEqn.H"
        }

        turbulence->correct();

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

    return(0);
}


// ************************************************************************* //
