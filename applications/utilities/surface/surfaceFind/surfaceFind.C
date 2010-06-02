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

Description
    Finds nearest triangle and vertex.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validOptions.insert("x", "X");
    argList::validOptions.insert("y", "Y");
    argList::validOptions.insert("z", "Z");

    argList::validArgs.append("surface file");

    argList args(argc, argv);

    point samplePt
    (
        args.optionRead<scalar>("x"),
        args.optionRead<scalar>("y"),
        args.optionRead<scalar>("z")
    );
    Info<< "Looking for nearest face/vertex to " << samplePt << endl;


    Info<< "Reading surf1 ..." << endl;
    triSurface surf1(args.additionalArgs()[0]);

    //
    // Nearest vertex
    //

    const pointField& localPoints = surf1.localPoints();

    label minIndex = -1;
    scalar minDist = GREAT;

    forAll(localPoints, pointI)
    {
        const scalar dist = mag(localPoints[pointI] - samplePt);
        if (dist < minDist)
        {
            minDist = dist;
            minIndex = pointI;
        }
    }

    Info<< "Nearest vertex:" << endl
        << "    index      :" << minIndex << " (in localPoints)" << endl
        << "    index      :" << surf1.meshPoints()[minIndex]
        << " (in points)" << endl
        << "    coordinates:" << localPoints[minIndex] << endl
        << endl;

    //
    // Nearest face
    //

    const pointField& points = surf1.points();

    minIndex = -1;
    minDist = GREAT;

    forAll(surf1, faceI)
    {
        const labelledTri& f = surf1[faceI];
        const point centre = f.centre(points);

        const scalar dist = mag(centre - samplePt);
        if (dist < minDist)
        {
            minDist = dist;
            minIndex = faceI;
        }
    }

    const labelledTri& f = surf1[minIndex];

    Info<< "Face with nearest centre:" << endl
        << "    index        :" << minIndex << endl
        << "    centre       :" << f.centre(points) << endl
        << "    face         :" << f << endl
        << "    vertex coords:" << points[f[0]] << " "
        << points[f[1]] << " " << points[f[2]] << endl
        << endl;


    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
