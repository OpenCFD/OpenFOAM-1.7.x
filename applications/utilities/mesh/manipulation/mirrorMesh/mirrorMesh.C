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

Description
    Mirrors a mesh around a given plane.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "mirrorFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    mirrorFvMesh mesh
    (
        IOobject
        (
            mirrorFvMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << "Writing mirrored mesh" << endl;
    mesh.mirrorMesh().write();

    Info << "End" << endl;

    return 0;
}


// ************************************************************************* //
