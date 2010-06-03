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
    Set normal consistent with respect to a user provided 'outside' point.
    If -inside the point is considered inside.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "orientedSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("Foam surface file");
    argList::validArgs.append("visiblePoint");
    argList::validArgs.append("output file");
    argList::validOptions.insert("inside", "");
    argList args(argc, argv);

    fileName surfFileName(args.additionalArgs()[0]);
    Info<< "Reading surface from " << surfFileName << endl;

    point visiblePoint(IStringStream(args.additionalArgs()[1])());
    Info<< "Visible point " << visiblePoint << endl;

    bool orientInside = args.optionFound("inside");

    if (orientInside)
    {
        Info<< "Orienting surface such that visiblePoint " << visiblePoint
            << " is inside" << endl;
    }
    else
    {
        Info<< "Orienting surface such that visiblePoint " << visiblePoint
            << " is outside" << endl;
    }

    fileName outFileName(args.additionalArgs()[2]);
    Info<< "Writing surface to " << outFileName << endl;


    // Load surface
    triSurface surf(surfFileName);

    //orientedSurface normalSurf(surf, visiblePoint, !orientInside);
    bool anyFlipped = orientedSurface::orient
    (
        surf,
        visiblePoint,
       !orientInside
    );

    if (anyFlipped)
    {
        Info<< "Flipped orientation of (part of) surface." << endl;
    }
    else
    {
        Info<< "Did not flip orientation of any triangle of surface." << endl;
    }

    Info<< "Writing new surface to " << outFileName << endl;

    surf.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
