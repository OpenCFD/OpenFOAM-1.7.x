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

\*---------------------------------------------------------------------------*/

//#include "octreeDataPointTreeLeaf.H"
#include "octreeDataPoint.H"
#include "treeLeaf.H"

// * * * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * //

template<>
Foam::label Foam::treeLeaf<Foam::octreeDataPoint>::find
(
    const octreeDataPoint& shapes,
    const point& sample
) const
{
    notImplemented
    (
        "Foam::treeLeaf<Foam::octreeDataPoint>::find("
        "const octreeDataPoint& shapes,"
        "const point& sample"
    );

    return false;
}


template<>
bool Foam::treeLeaf<Foam::octreeDataPoint>::findNearest
(
    const octreeDataPoint& shapes,
    const point& sample,
    treeBoundBox& tightest,
    label& tightestI,
    scalar& tightestDist
) const
{
    // Some aliases
    const pointField& points = shapes.points();
    point& tMin = tightest.min();
    point& tMax = tightest.max();

    scalar minDist2 = sqr(tightestDist);

    label minIndex = -1;
    forAll(indices_, i)
    {
        label pointi = indices_[i];
        scalar dist = magSqr(points[pointi] - sample);

        if (dist < minDist2)
        {
            minDist2 = dist;
            minIndex = pointi;
        }
    }

    if (minIndex != -1)
    {
        tightestDist = sqrt(minDist2);
    
        // New nearer. Update 'tightest' bounding box
        tMin.x() = sample.x() - tightestDist;
        tMin.y() = sample.y() - tightestDist;
        tMin.z() = sample.z() - tightestDist;

        tMax.x() = sample.x() + tightestDist;
        tMax.y() = sample.y() + tightestDist;
        tMax.z() = sample.z() + tightestDist;

        tightestI = minIndex;

        return true;
    }
    else
    {
        // New no nearer so nothing changed
        return false;
    }
}


// ************************************************************************* //
