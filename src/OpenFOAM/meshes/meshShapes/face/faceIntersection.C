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

Description
    Return intersection of a line with the face

\*---------------------------------------------------------------------------*/

#include "face.H"
#include "pointHit.H"
#include "triPointRef.H"
#include "line.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return potential intersection with face with a ray starting
// at p, direction n (does not need to be normalized)
// Does face-center decomposition and returns triangle intersection
// point closest to p.

// In case of miss the point is the nearest point intersection of the
// face plane and the ray and the distance is the distance between the
// intersection point and the nearest point on the face

Foam::pointHit Foam::face::ray
(
    const point& p,
    const vector& n,
    const pointField& meshPoints,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    point ctr = Foam::average(points(meshPoints));

    scalar nearestHitDist = GREAT;

    scalar nearestMissDist = GREAT;
    bool eligible = false;

    // Initialize to miss, distance = GREAT
    pointHit nearest(p);

    const labelList& f = *this;

    label nPoints = size();

    point nextPoint = ctr;

    for (label pI = 0; pI < nPoints; pI++)
    {
        nextPoint = meshPoints[f[fcIndex(pI)]];

        // Note: for best accuracy, centre point always comes last
        //
        pointHit curHit = triPointRef
        (
            meshPoints[f[pI]],
            nextPoint,
            ctr
        ).ray(p, n, alg, dir);

        if (curHit.hit())
        {
            if (Foam::mag(curHit.distance()) < Foam::mag(nearestHitDist))
            {
                nearestHitDist = curHit.distance();
                nearest.setHit();
                nearest.setPoint(curHit.hitPoint());
            }
        }
        else if (!nearest.hit())
        {
            // Miss and no hit yet. Update miss statistics.
            if (curHit.eligibleMiss())
            {
                eligible = true;

                // Miss distance is the distance between the plane intersection
                // point and the nearest point of the triangle
                scalar missDist =
                    Foam::mag
                    (
                        p + curHit.distance()*n
                      - curHit.missPoint()
                    );

                if (missDist < nearestMissDist)
                {
                    nearestMissDist = missDist;
                    nearest.setDistance(curHit.distance());
                    nearest.setPoint(curHit.missPoint());
                }
            }
        }
    }

    if (nearest.hit())
    {
        nearest.setDistance(nearestHitDist);
    }
    else
    {
        // Haven't hit a single face triangle
        nearest.setMiss(eligible);
    }

    return nearest;
}


Foam::pointHit Foam::face::intersection
(
    const point& p,
    const vector& q,
    const point& ctr,
    const pointField& meshPoints,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    scalar nearestHitDist = VGREAT;

    // Initialize to miss, distance = GREAT
    pointHit nearest(p);

    const labelList& f = *this;

    forAll(f, pI)
    {
        // Note: for best accuracy, centre point always comes last
        pointHit curHit = triPointRef
        (
            meshPoints[f[pI]],
            meshPoints[f[fcIndex(pI)]],
            ctr
        ).intersection(p, q, alg, tol);

        if (curHit.hit())
        {
            if (Foam::mag(curHit.distance()) < nearestHitDist)
            {
                nearestHitDist = Foam::mag(curHit.distance());
                nearest.setHit();
                nearest.setPoint(curHit.hitPoint());
            }
        }
    }

    if (nearest.hit())
    {
        nearest.setDistance(nearestHitDist);
    }

    return nearest;
}


Foam::pointHit Foam::face::nearestPoint
(
    const point& p,
    const pointField& meshPoints
) const
{
    const face& f = *this;
    point ctr = centre(meshPoints);

    // Initialize to miss, distance=GREAT
    pointHit nearest(p);

    label nPoints = f.size();

    point nextPoint = ctr;

    for (label pI = 0; pI < nPoints; pI++)
    {
        nextPoint = meshPoints[f[fcIndex(pI)]];

        // Note: for best accuracy, centre point always comes last
        //
        triPointRef tri
        (
            meshPoints[f[pI]],
            nextPoint,
            ctr
        );

        pointHit curHit = tri.nearestPoint(p);

        if (Foam::mag(curHit.distance()) < Foam::mag(nearest.distance()))
        {
            nearest.setDistance(curHit.distance());

            if (curHit.hit())
            {
                nearest.setHit();
                nearest.setPoint(curHit.hitPoint());
            }
            else
            {
                // In nearest point, miss is always eligible
                nearest.setMiss(true);
                nearest.setPoint(curHit.missPoint());
            }
        }
    }

    return nearest;
}


// ************************************************************************* //
