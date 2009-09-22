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

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool triSurface::stitchTriangles
(
    const pointField& rawPoints,
    const scalar tol,
    bool verbose
)
{
    // Merge points
    labelList pointMap;
    pointField newPoints;
    bool hasMerged = mergePoints(rawPoints, tol, verbose, pointMap, newPoints);

    pointField& ps = storedPoints();

    // Set the coordinates to the merged ones
    ps = newPoints;

    if (hasMerged)
    {
        if (verbose)
        {
            Pout<< "stitchTriangles : Merged from " << rawPoints.size()
                << " points down to " << newPoints.size() << endl;
        }

        // Reset the triangle point labels to the unique points array
        label newTriangleI = 0;
        forAll(*this, i)
        {
            label newA = pointMap[operator[](i)[0]];
            label newB = pointMap[operator[](i)[1]];
            label newC = pointMap[operator[](i)[2]];

            if ((newA != newB) && (newA != newC) && (newB != newC))
            {
                operator[](newTriangleI)[0] = newA;
                operator[](newTriangleI)[1] = newB;
                operator[](newTriangleI)[2] = newC;
                operator[](newTriangleI).region() = operator[](i).region();
                newTriangleI++;
            }
            else if (verbose)
            {
                Pout<< "stitchTriangles : "
                    << "Removing triangle " << i << " with non-unique vertices."
                    << endl
                    << "    vertices   :" << newA << ' ' << newB << ' ' << newC
                    << endl
                    << "    coordinates:" << ps[newA] << ' ' << ps[newB]
                    << ' ' << ps[newC]  << endl;
            }
        }

        if (newTriangleI != size())
        {
            if (verbose)
            {
                Pout<< "stitchTriangles : "
                    << "Removed " << size() - newTriangleI
                    << " triangles" << endl;
            }
            setSize(newTriangleI);
        }
    }

    return hasMerged;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
