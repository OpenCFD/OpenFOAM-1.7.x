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
    Example of simple laplacian smoother

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"
#include "boundBox.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.clear();
    argList::validArgs.clear();
    argList::validArgs.append("surface file");
    argList::validArgs.append("underrelax factor (0..1)");
    argList::validArgs.append("iterations");
    argList::validArgs.append("output file");
    argList args(argc, argv);

    fileName surfFileName(args.additionalArgs()[0]);
    scalar relax(readScalar(IStringStream(args.additionalArgs()[1])()));
    if ((relax <= 0) || (relax > 1))
    {
        FatalErrorIn(args.executable()) << "Illegal relaxation factor "
            << relax << endl
            << "0: no change   1: move vertices to average of neighbours"
            << exit(FatalError);
    }
    label iters(readLabel(IStringStream(args.additionalArgs()[2])()));
    fileName outFileName(args.additionalArgs()[3]);

    Info<< "Relax:" << relax << endl;
    Info<< "Iters:" << iters << endl;


    Info<< "Reading surface from " << surfFileName << " ..." << endl;

    triSurface surf1(surfFileName);

    Info<< "Triangles    : " << surf1.size() << endl;
    Info<< "Vertices     : " << surf1.nPoints() << endl;
    Info<< "Bounding Box : " << boundBox(surf1.localPoints()) << endl;

    pointField newPoints(surf1.localPoints());

    const labelListList& pointEdges = surf1.pointEdges();


    for(label iter = 0; iter < iters; iter++)
    {
        forAll(pointEdges, vertI)
        {
            vector avgPos(vector::zero);

            const labelList& myEdges = pointEdges[vertI];

            forAll(myEdges, myEdgeI)
            {
                const edge& e = surf1.edges()[myEdges[myEdgeI]];

                label otherVertI = e.otherVertex(vertI);

                avgPos += surf1.localPoints()[otherVertI];
            }
            avgPos /= myEdges.size();

            newPoints[vertI] = (1-relax)*newPoints[vertI] + relax*avgPos;
        }
    }

    triSurface surf2
    (
        surf1.localFaces(),
        surf1.patches(),
        newPoints
    );

    Info<< "Writing surface to " << outFileName << " ..." << endl;

    surf2.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
