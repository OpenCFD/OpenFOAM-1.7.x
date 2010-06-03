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

\*---------------------------------------------------------------------------*/

#include "mergePoints.H"
#include "SortableList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mergePoints
(
    const UList<point>& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    List<point>& newPoints,
    const point& origin
)
{
    point compareOrigin = origin;

    if (origin == point(VGREAT, VGREAT, VGREAT))
    {
        if (points.size())
        {
            compareOrigin = sum(points)/points.size();
        }
    }

    // Create a old to new point mapping array
    pointMap.setSize(points.size());
    pointMap = -1;

    // Storage for merged points
    newPoints.setSize(points.size());

    if (points.empty())
    {
        return false;
    }


    const scalar mergeTolSqr = sqr(mergeTol);

    // Sort points by magSqr
    SortableList<scalar> sortedMagSqr(magSqr(points - compareOrigin));

    bool hasMerged = false;

    label newPointI = 0;


    // Handle 0th point separately (is always unique)
    label pointI = sortedMagSqr.indices()[0];
    pointMap[pointI] = newPointI;
    newPoints[newPointI++] = points[pointI];


    for (label sortI = 1; sortI < sortedMagSqr.size(); sortI++)
    {
        // Get original point index
        label pointI = sortedMagSqr.indices()[sortI];

        // Compare to previous points to find equal one.
        label equalPointI = -1;

        for
        (
            label prevSortI = sortI - 1;
            prevSortI >= 0
         && mag
            (
                sortedMagSqr[prevSortI]
             -  sortedMagSqr[sortI]
            ) <= mergeTolSqr;
            prevSortI--
        )
        {
            label prevPointI = sortedMagSqr.indices()[prevSortI];

            if (magSqr(points[pointI] - points[prevPointI]) <= mergeTolSqr)
            {
                // Found match.
                equalPointI = prevPointI;

                break;
            }
        }


        if (equalPointI != -1)
        {
            // Same coordinate as equalPointI. Map to same new point.
            pointMap[pointI] = pointMap[equalPointI];

            hasMerged = true;

            if (verbose)
            {
                Pout<< "Foam::mergePoints : Merging points "
                    << pointI << " and " << equalPointI
                    << " with coordinates:" << points[pointI]
                    << " and " << points[equalPointI]
                    << endl;
            }
        }
        else
        {
            // Differs. Store new point.
            pointMap[pointI] = newPointI;
            newPoints[newPointI++] = points[pointI];
        }
    }

    newPoints.setSize(newPointI);

    return hasMerged;
}


// ************************************************************************* //
