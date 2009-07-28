/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    dsmcFoam

Description
    Direct simulation Monte Carlo (DSMC) solver for 3D, transient, multi-
    species flows

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Carry out dsmcCloud timestep

        dsmc.evolve();

        // Retrieve flow field data from dsmcCloud

        rhoN = dsmc.rhoN();
        rhoN.correctBoundaryConditions();

        rhoM = dsmc.rhoM();
        rhoM.correctBoundaryConditions();

        dsmcRhoN = dsmc.dsmcRhoN();
        dsmcRhoN.correctBoundaryConditions();

        momentum = dsmc.momentum();
        momentum.correctBoundaryConditions();

        linearKE = dsmc.linearKE();
        linearKE.correctBoundaryConditions();

        internalE = dsmc.internalE();
        internalE.correctBoundaryConditions();

        iDof = dsmc.iDof();
        iDof.correctBoundaryConditions();

        // Retrieve surface field data from dsmcCloud

        q = dsmc.q();

        fD = dsmc.fD();

        // Print status of dsmcCloud

        dsmc.info();

        runTime.write();

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
