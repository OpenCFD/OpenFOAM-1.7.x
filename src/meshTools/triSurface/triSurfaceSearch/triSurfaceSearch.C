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

#include "triSurfaceSearch.H"
#include "indexedOctree.H"
#include "boolList.H"
#include "treeDataTriSurface.H"
#include "triSurface.H"
#include "line.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const point triSurfaceSearch::greatPoint(GREAT, GREAT, GREAT);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface. Holds reference!
triSurfaceSearch::triSurfaceSearch(const triSurface& surface)
:
    surface_(surface),
    treePtr_(NULL)
{
    // Random number generator. Bit dodgy since not exactly random ;-)
    Random rndGen(65431);

    // Slightly extended bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    treeBoundBox treeBb
    (
        treeBoundBox(surface_.points(), surface_.meshPoints()).extend
        (
            rndGen,
            1E-4
        )
    );
    treeBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    treeBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    treePtr_.reset
    (
        new indexedOctree<treeDataTriSurface>
        (
            treeDataTriSurface(surface_),
            treeBb,
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Determine inside/outside for samples
boolList triSurfaceSearch::calcInside(const pointField& samples) const
{
    boolList inside(samples.size());

    forAll(samples, sampleI)
    {
        const point& sample = samples[sampleI];

        if (!tree().bb().contains(sample))
        {
            inside[sampleI] = false;
        }
        else if
        (
            tree().getVolumeType(sample)
         == indexedOctree<treeDataTriSurface>::INSIDE
        )
        {
            inside[sampleI] = true;
        }
        else
        {
            inside[sampleI] = false;
        }
    }
    return inside;
}


labelList triSurfaceSearch::calcNearestTri
(
    const pointField& samples,
    const vector& span
) const
{
    labelList nearest(samples.size());

    const scalar nearestDistSqr = 0.25*magSqr(span);

    pointIndexHit hitInfo;

    forAll(samples, sampleI)
    {
        hitInfo = tree().findNearest(samples[sampleI], nearestDistSqr);

        if (hitInfo.hit())
        {
            nearest[sampleI] = hitInfo.index();
        }
        else
        {
            nearest[sampleI] = -1;
        }
    }

    return nearest;
}


// Nearest point on surface
tmp<pointField> triSurfaceSearch::calcNearest
(
    const pointField& samples,
    const vector& span
) const
{
    const scalar nearestDistSqr = 0.25*magSqr(span);

    tmp<pointField> tnearest(new pointField(samples.size()));
    pointField& nearest = tnearest();

    pointIndexHit hitInfo;

    forAll(samples, sampleI)
    {
        hitInfo = tree().findNearest(samples[sampleI], nearestDistSqr);

        if (hitInfo.hit())
        {
            nearest[sampleI] = hitInfo.hitPoint();
        }
        else
        {
            nearest[sampleI] = greatPoint;
        }
    }

    return tnearest;
}


pointIndexHit triSurfaceSearch::nearest(const point& pt, const vector& span)
 const
{
    const scalar nearestDistSqr = 0.25*magSqr(span);

    return tree().findNearest(pt, nearestDistSqr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
