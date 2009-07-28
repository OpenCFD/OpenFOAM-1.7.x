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

Description
    Calculate distance to wall.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallDistData.H"
#include "wallPointYPlus.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info<< "Mesh read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;


    Info<< "Time now = " << runTime.timeName() << endl;

    // wall distance and yStar

    volScalarField yStar
    (
        IOobject
        (
            "yStar",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("yStar", dimless, 1.0)
    );

    // Fill wall patches of yStar with some value.
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (typeid(patch) == typeid(wallFvPatch))
        {
            fvPatchScalarField& wallData = yStar.boundaryField()[patchI];

            forAll(patch, patchFaceI)
            {
// Hack. Just some value.
                wallData[patchFaceI] = 1/2500.0;
            }
        }
    }


    // Do distance calculation, transporting values of yStar
    wallPointYPlus::yPlusCutOff = 200;
    wallDistData<wallPointYPlus> y(mesh, yStar, true);

    if (y.nUnset() != 0)
    {
        WarningIn(args.executable())
            << "There are " << y.nUnset()
            << " remaining unset cells and/or boundary values" << endl;
    }


    y.write();

    y.data().write();

    volScalarField yPlus
    (
        IOobject
        (
            "yPlus",
            mesh.time().timeName(),
            mesh
        ),
        y/y.data()
    );

    yPlus.write();

    Info<< "End\n" << endl;

    return 0;



}


// ************************************************************************* //
