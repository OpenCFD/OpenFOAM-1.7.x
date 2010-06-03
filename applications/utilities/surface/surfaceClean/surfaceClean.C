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
    surfaceClean

Description
    - collapses small edges, removing triangles.
    - converts sliver triangles into split edges by projecting point onto
      base of triangle.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"

#include "collapseBase.H"
#include "collapseEdge.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList::noParallel();
    argList::validArgs.append("surface file");
    argList::validArgs.append("min length");
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName inFileName(args.additionalArgs()[0]);
    scalar minLen(readScalar(IStringStream(args.additionalArgs()[1])()));
    fileName outFileName(args.additionalArgs()[2]);

    Pout<< "Reading surface " << inFileName << nl
        << "Collapsing all triangles with edges or heights < " << minLen << nl
        << "Writing result to " << outFileName << nl << endl;


    Pout<< "Reading surface from " << inFileName << " ..." << nl << endl;
    triSurface surf(inFileName);
    surf.writeStats(Pout);


    Pout<< "Collapsing triangles to edges ..." << nl << endl;

    while (true)
    {
        label nEdgeCollapse = collapseEdge(surf, minLen);

        if (nEdgeCollapse == 0)
        {
            break;
        }
    }
    while (true)
    {
        label nSplitEdge = collapseBase(surf, minLen);

        if (nSplitEdge == 0)
        {
            break;
        }
    }

    Pout<< nl << "Resulting surface:" << endl;
    surf.writeStats(Pout);
    Pout<< nl;

    Pout<< "Writing refined surface to " << outFileName << " ..." << endl;
    surf.write(outFileName);
    Pout<< nl;

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
