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

#include "searchableCylinder.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableCylinder, 0);
addToRunTimeSelectionTable(searchableSurface, searchableCylinder, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableCylinder::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    vector v(sample - point1_);

    // Decompose sample-point1 into normal and parallel component
    scalar parallel = (v & unitDir_);

    // Remove the parallel component
    v -= parallel*unitDir_;
    scalar magV = mag(v);

    if (parallel <= 0)
    {
        // nearest is at point1 end cap. Clip to radius.
        if (magV < ROOTVSMALL)
        {
            info.setPoint(point1_);
        }
        else
        {
            //info.setPoint(point1_ + min(magV, radius_)*v/magV);
            info.setPoint(point1_ + radius_*(v/magV));
        }
    }
    else if (parallel >= magDir_)
    {
        // nearest is at point2
        if (magV < ROOTVSMALL)
        {
            info.setPoint(point2_);
        }
        else
        {
            info.setPoint(point2_ + min(magV, radius_)*v/magV);
        }
    }
    else
    {
        if (magV < ROOTVSMALL)
        {
            info.setPoint(point1_ + parallel*unitDir_);
        }
        else
        {
            // Project onto radius
            info.setPoint(point1_ + parallel*unitDir_ + radius_*v/magV);
        }
    }

    if (magSqr(sample - info.rawPoint()) < nearestDistSqr)
    {
        info.setHit();
        info.setIndex(0);
    }

    return info;
}


Foam::scalar Foam::searchableCylinder::radius2(const point& pt) const
{
    const vector x = (pt-point1_) ^ unitDir_;
    return x&x;
}


// From http://www.gamedev.net/community/forums/topic.asp?topic_id=467789 -
// intersection of cylinder with ray
void Foam::searchableCylinder::findLineAll
(
    const point& start,
    const point& end,
    pointIndexHit& near,
    pointIndexHit& far
) const
{
    near.setMiss();
    far.setMiss();

    vector point1Start(start-point1_);
    vector point2Start(start-point2_);
    vector point1End(end-point1_);

    // Quick rejection of complete vector outside endcaps
    scalar s1 = point1Start&unitDir_;
    scalar s2 = point1End&unitDir_;

    if ((s1 < 0 && s2 < 0) || (s1 > magDir_ && s2 > magDir_))
    {
        return;
    }

    // Line as P = start+t*V  where V is unit vector and t=[0..mag(end-start)]
    vector V(end-start);
    scalar magV = mag(V);
    if (magV < ROOTVSMALL)
    {
        return;
    }
    V /= magV;


    // We now get the nearest intersections to start. This can either be
    // the intersection with the end plane or with the cylinder side.

    // Get the two points (expressed in t) on the end planes. This is to
    // clip any cylinder intersection against.
    scalar tPoint1;
    scalar tPoint2;

    // Maintain the two intersections with the endcaps
    scalar tNear = VGREAT;
    scalar tFar = VGREAT;

    {
        scalar s = (V&unitDir_);
        if (mag(s) > VSMALL)
        {
            tPoint1 = -s1/s;
            tPoint2 = -(point2Start&unitDir_)/s;
            if (tPoint2 < tPoint1)
            {
                Swap(tPoint1, tPoint2);
            }
            if (tPoint1 > magV || tPoint2 < 0)
            {
                return;
            }

            // See if the points on the endcaps are actually inside the cylinder
            if (tPoint1 >= 0 && tPoint1 <= magV)
            {
                if (radius2(start+tPoint1*V) <= sqr(radius_))
                {
                    tNear = tPoint1;
                }
            }
            if (tPoint2 >= 0 && tPoint2 <= magV)
            {
                if (radius2(start+tPoint2*V) <= sqr(radius_))
                {
                    // Check if already have a near hit from point1
                    if (tNear <= 1)
                    {
                        tFar = tPoint2;
                    }
                    else
                    {
                        tNear = tPoint2;
                    }
                }
            }
        }
        else
        {
            // Vector perpendicular to cylinder. Check for outside already done
            // above so just set tpoint to allow all.
            tPoint1 = -VGREAT;
            tPoint2 = VGREAT;
        }
    }


    const vector x = point1Start ^ unitDir_;
    const vector y = V ^ unitDir_;
    const scalar d = sqr(radius_);

    // Second order equation of the form a*t^2 + b*t + c
    const scalar a = (y&y);
    const scalar b = 2*(x&y);
    const scalar c = (x&x)-d;

    const scalar disc = b*b-4*a*c;

    scalar t1 = -VGREAT;
    scalar t2 = VGREAT;

    if (disc < 0)
    {
        // Fully outside
        return;
    }
    else if (disc < ROOTVSMALL)
    {
        // Single solution
        if (mag(a) > ROOTVSMALL)
        {
            t1 = -b/(2*a);

            //Pout<< "single solution t:" << t1
            //    << " for start:" << start << " end:" << end
            //    << " c:" << c << endl;

            if (t1 >= 0 && t1 <= magV && t1 >= tPoint1 && t1 <= tPoint2)
            {
                // valid. Insert sorted.
                if (t1 < tNear)
                {
                    tFar = tNear;
                    tNear = t1;
                }
                else if (t1 < tFar)
                {
                    tFar = t1;
                }
            }
            else
            {
                return;
            }
        }
        else
        {
            // Aligned with axis. Check if outside radius
            //Pout<< "small discriminant:" << disc
            //    << " for start:" << start << " end:" << end
            //    << " magV:" << magV
            //    << " c:" << c << endl;
            if (c > 0)
            {
                return;
            }
        }
    }
    else
    {
        if (mag(a) > ROOTVSMALL)
        {
            scalar sqrtDisc = sqrt(disc);

            t1 = (-b - sqrtDisc)/(2*a);
            t2 = (-b + sqrtDisc)/(2*a);
            if (t2 < t1)
            {
                Swap(t1, t2);
            }

            if (t1 >= 0 && t1 <= magV && t1 >= tPoint1 && t1 <= tPoint2)
            {
                // valid. Insert sorted.
                if (t1 < tNear)
                {
                    tFar = tNear;
                    tNear = t1;
                }
                else if (t1 < tFar)
                {
                    tFar = t1;
                }
            }
            if (t2 >= 0 && t2 <= magV && t2 >= tPoint1 && t2 <= tPoint2)
            {
                // valid. Insert sorted.
                if (t2 < tNear)
                {
                    tFar = tNear;
                    tNear = t2;
                }
                else if (t2 < tFar)
                {
                    tFar = t2;
                }
            }
            //Pout<< "two solutions t1:" << t1 << " t2:" << t2
            //    << " for start:" << start << " end:" << end
            //    << " magV:" << magV
            //    << " c:" << c << endl;
        }
        else
        {
            // Aligned with axis. Check if outside radius
            //Pout<< "large discriminant:" << disc
            //    << " small a:" << a
            //    << " for start:" << start << " end:" << end
            //    << " magV:" << magV
            //    << " c:" << c << endl;
            if (c > 0)
            {
                return;
            }
        }
    }

    // Check tNear, tFar
    if (tNear >= 0 && tNear <= magV)
    {
        near.setPoint(start+tNear*V);
        near.setHit();
        near.setIndex(0);

        if (tFar <= magV)
        {
            far.setPoint(start+tFar*V);
            far.setHit();
            far.setIndex(0);
        }
    }
    else if (tFar >= 0 && tFar <= magV)
    {
        near.setPoint(start+tFar*V);
        near.setHit();
        near.setIndex(0);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableCylinder::searchableCylinder
(
    const IOobject& io,
    const point& point1,
    const point& point2,
    const scalar radius
)
:
    searchableSurface(io),
    point1_(point1),
    point2_(point2),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_),
    radius_(radius)
{}


Foam::searchableCylinder::searchableCylinder
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    point1_(dict.lookup("point1")),
    point2_(dict.lookup("point2")),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableCylinder::~searchableCylinder()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableCylinder::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


void Foam::searchableCylinder::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableCylinder::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected take second one.
        pointIndexHit b;
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableCylinder::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        // Discard far intersection
        pointIndexHit b;
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableCylinder::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit> >& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;
        findLineAll(start[i], end[i], near, far);

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].setSize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].setSize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].setSize(1);
                info[i][0] = far;
            }
            else
            {
                info[i].clear();
            }
        }
    }
}


void Foam::searchableCylinder::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableCylinder::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = vector::zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            vector v(info[i].hitPoint() - point1_);

            // Decompose sample-point1 into normal and parallel component
            scalar parallel = v & unitDir_;

            if (parallel < 0)
            {
                normal[i] = -unitDir_;
            }
            else if (parallel > magDir_)
            {
                normal[i] = -unitDir_;
            }
            else
            {
                // Remove the parallel component
                v -= parallel*unitDir_;
                normal[i] = v/mag(v);
            }
        }
    }
}


void Foam::searchableCylinder::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());
    volType = INSIDE;

    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        vector v(pt - point1_);

        // Decompose sample-point1 into normal and parallel component
        scalar parallel = v & unitDir_;

        if (parallel < 0)
        {
            // left of point1 endcap
            volType[pointI] = OUTSIDE;
        }
        else if (parallel > magDir_)
        {
            // right of point2 endcap
            volType[pointI] = OUTSIDE;
        }
        else
        {
            // Remove the parallel component
            v -= parallel*unitDir_;

            if (mag(v) > radius_)
            {
                volType[pointI] = OUTSIDE;
            }
            else
            {
                volType[pointI] = INSIDE;
            }
        }
    }
}


// ************************************************************************* //
