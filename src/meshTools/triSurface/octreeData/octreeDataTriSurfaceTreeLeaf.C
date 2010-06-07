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

\*---------------------------------------------------------------------------*/

#include "octreeDataTriSurfaceTreeLeaf.H"
#include "octreeDataTriSurface.H"

// * * * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * //

template<>
bool Foam::treeLeaf<Foam::octreeDataTriSurface>::findNearest
(
    const octreeDataTriSurface& shapes,
    const point& sample,
    treeBoundBox& tightest,
    label& tightestI,
    scalar& tightestDist
) const
{
    // Some aliases
    const treeBoundBoxList& allBb = shapes.allBb();
    point& min = tightest.min();
    point& max = tightest.max();

    point nearest;

    bool changed = false;
    forAll(indices_, i)
    {
        label faceI = indices_[i];

        // Quick rejection test.
        if (tightest.overlaps(allBb[faceI]))
        {
            // Full calculation
            scalar dist = shapes.calcNearest(faceI, sample, nearest);

            if (dist < tightestDist)
            {
                // Update bb (centered around sample, span is dist)
                min.x() = sample.x() - dist;
                min.y() = sample.y() - dist;
                min.z() = sample.z() - dist;

                max.x() = sample.x() + dist;
                max.y() = sample.y() + dist;
                max.z() = sample.z() + dist;

                tightestI = faceI;
                tightestDist = dist;

                changed = true;
            }
        }
    }
    return changed;
}


// ************************************************************************* //
