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
    engineCompRatioFoam

Description
    Calculate the geometric compression ratio.
    Note that if you have valves and/or extra volumes it will not work,
    since it calculates the volume at BDC and TCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "engineTime.H"
#include "engineMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createEngineTime.H"
    #include "createEngineMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar eps = 1.0e-10;
    scalar fullCycle = 360.0;

    scalar ca0 = -180.0;
    scalar ca1 = 0.0;

    while (runTime.theta() > ca0)
    {
        ca0 += fullCycle;
        ca1 += fullCycle;
    }

    while (mag(runTime.theta() - ca0) > eps)
    {
        scalar t0 = runTime.userTimeToTime(ca0 - runTime.theta());
        runTime.setDeltaT(t0);
        runTime++;
        Info<< "CA = " << runTime.theta() << endl;
        mesh.move();
    }

    scalar Vmax = sum(mesh.V().field());

    while (mag(runTime.theta()-ca1) > eps)
    {
        scalar t1 = runTime.userTimeToTime(ca1-runTime.theta());
        runTime.setDeltaT(t1);
        runTime++;
        Info<< "CA = " << runTime.theta() << endl;
        mesh.move();
    }

    scalar Vmin = sum(mesh.V().field());

    Info<< "\nVmax = " << Vmax;
    Info<< ", Vmin = " << Vmin << endl;
    Info<< "Vmax/Vmin = " << Vmax/Vmin << endl;
    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
