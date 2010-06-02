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

#include "octreeDataTriSurface.H"

#include "labelList.H"
#include "treeBoundBox.H"
#include "faceList.H"
#include "triPointRef.H"
#include "octree.H"
#include "triSurfaceTools.H"
#include "triangleFuncs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::octreeDataTriSurface, 0);

Foam::scalar Foam::octreeDataTriSurface::tol(1E-6);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Fast distance to triangle calculation. From
// "Distance Between Point and Trangle in 3D"
// David Eberly, Magic Software Inc. Aug. 2003.
// Works on function Q giving distance to point and tries to minimize this.
void Foam::octreeDataTriSurface::nearestCoords
(
    const point& base,
    const point& E0,
    const point& E1,
    const scalar a,
    const scalar b,
    const scalar c,
    const point& P,
    scalar& s,
    scalar& t
)
{
    // distance vector
    const vector D(base - P);

    // Precalculate distance factors.
    const scalar d = E0 & D;
    const scalar e = E1 & D;

    // Do classification
    const scalar det = a*c - b*b;

    s = b*e - c*d;
    t = b*d - a*e;

    if (s+t < det)
    {
        if (s < 0)
        {
            if (t < 0)
            {
                //region 4
                if (e > 0)
                {
                    //min on edge t = 0
                    t = 0;
                    s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
                }
                else
                {
                    //min on edge s=0
                    s = 0;
                    t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
                }
            }
            else
            {
                //region 3. Min on edge s = 0
                s = 0;
                t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
            }
        }
        else if (t < 0)
        {
            //region 5
            t = 0;
            s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
        }
        else
        {
            //region 0
            const scalar invDet = 1/det;
            s *= invDet;
            t *= invDet;
        }
    }
    else
    {
        if (s < 0)
        {
            //region 2
            const scalar tmp0 = b + d;
            const scalar tmp1 = c + e;
            if (tmp1 > tmp0)
            {
                //min on edge s+t=1
                const scalar numer = tmp1 - tmp0;
                const scalar denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
                t = 1 - s;
            }
            else
            {
                //min on edge s=0
                s = 0;
                t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : - e/c));
            }
        }
        else if (t < 0)
        {
            //region 6
            const scalar tmp0 = b + d;
            const scalar tmp1 = c + e;
            if (tmp1 > tmp0)
            {
                //min on edge s+t=1
                const scalar numer = tmp1 - tmp0;
                const scalar denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
                t = 1 - s;
            }
            else
            {
                //min on edge t=0
                t = 0;
                s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : - d/a));
            }
        }
        else
        {
            //region 1
            const scalar numer = c+e-(b+d);
            if (numer <= 0)
            {
                s = 0;
            }
            else
            {
                const scalar denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
            }
        }
        t = 1 - s;
    }


    // Calculate distance.
    // Note: abs should not be needed but truncation error causes problems
    // with points very close to one of the triangle vertices.
    // (seen up to -9e-15). Alternatively add some small value.

    // const scalar f = D & D;
    // return a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f + SMALL;
    // return Foam::mag(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f);
}


Foam::point Foam::octreeDataTriSurface::nearestPoint
(
    const label index,
    const point& p
) const
{
    scalar s;
    scalar t;

    nearestCoords
    (
        base_[index],
        E0_[index],
        E1_[index],
        a_[index],
        b_[index],
        c_[index],
        p,
        s,
        t
    );

    return base_[index] + s*E0_[index] + t*E1_[index];
}


// Helper function to calculate tight fitting bounding boxes.
Foam::treeBoundBoxList Foam::octreeDataTriSurface::calcBb
(
    const triSurface& surf
)
{
    treeBoundBoxList allBb(surf.size(), treeBoundBox::invertedBox);

    const labelListList& pointFcs = surf.pointFaces();
    const pointField& localPts = surf.localPoints();

    forAll(pointFcs, pointI)
    {
        const labelList& myFaces = pointFcs[pointI];
        const point& vertCoord = localPts[pointI];

        forAll(myFaces, myFaceI)
        {
            // Update bb
            label faceI = myFaces[myFaceI];

            treeBoundBox& bb = allBb[faceI];

            bb.min() = min(bb.min(), vertCoord);
            bb.max() = max(bb.max(), vertCoord);
        }
    }

    return allBb;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::octreeDataTriSurface::octreeDataTriSurface(const triSurface& surface)
:
    surface_(surface),
    allBb_(calcBb(surface_)),
    base_(surface_.size()),
    E0_(surface_.size()),
    E1_(surface_.size()),
    a_(surface_.size()),
    b_(surface_.size()),
    c_(surface_.size())
{
    // Precalculate factors for distance calculation
    const pointField& points = surface_.points();

    forAll(surface_, faceI)
    {
        const labelledTri& f = surface_[faceI];

        // Calculate base and spanning vectors of triangles
        base_[faceI] = points[f[1]];
        E0_[faceI] = points[f[0]] - points[f[1]];
        E1_[faceI] = points[f[2]] - points[f[1]];

        a_[faceI] = E0_[faceI] & E0_[faceI];
        b_[faceI] = E0_[faceI] & E1_[faceI];
        c_[faceI] = E1_[faceI] & E1_[faceI];
    }
}


// Construct from components
Foam::octreeDataTriSurface::octreeDataTriSurface
(
    const triSurface& surface,
    const treeBoundBoxList& allBb
)
:
    surface_(surface),
    allBb_(allBb),
    base_(surface_.size()),
    E0_(surface_.size()),
    E1_(surface_.size()),
    a_(surface_.size()),
    b_(surface_.size()),
    c_(surface_.size())
{
    const pointField& points = surface_.points();

    forAll(surface_, faceI)
    {
        const labelledTri& f = surface_[faceI];

        // Calculate base and spanning vectors of triangles
        base_[faceI] = points[f[1]];
        E0_[faceI] = points[f[0]] - points[f[1]];
        E1_[faceI] = points[f[2]] - points[f[1]];

        a_[faceI] = E0_[faceI] & E0_[faceI];
        b_[faceI] = E0_[faceI] & E1_[faceI];
        c_[faceI] = E1_[faceI] & E1_[faceI];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::octreeDataTriSurface::getSampleType
(
    const octree<octreeDataTriSurface>& oc,
    const point& sample
) const

{
    treeBoundBox tightest(treeBoundBox::greatBox);
    scalar tightestDist(treeBoundBox::great);

    // Find nearest face to sample
    label faceI = oc.findNearest(sample, tightest, tightestDist);

    if (debug & 2)
    {
        Pout<< "getSampleType : sample:" << sample
            << " nearest face:" << faceI;
    }

    if (faceI == -1)
    {
        FatalErrorIn
        (
            "octreeDataTriSurface::getSampleType"
            "(octree<octreeDataTriSurface>&, const point&)"
        )   << "Could not find " << sample << " in octree."
            << abort(FatalError);
    }

    const pointField& pts = surface_.points();
    const labelledTri& f = surface_[faceI];

    pointHit curHit = triPointRef
    (
        pts[f[0]],
        pts[f[1]],
        pts[f[2]]
    ).nearestPoint(sample);

    // Get normal according to position on face. On point -> pointNormal,
    // on edge-> edge normal, face normal on interior.
    vector n
    (
        triSurfaceTools::surfaceNormal
        (
            surface_,
            faceI,
            curHit.rawPoint()
        )
    );

    return
        octree<octreeDataTriSurface>::getVolType(n, sample - curHit.rawPoint());
}


// Check if any point on triangle is inside cubeBb.
bool Foam::octreeDataTriSurface::overlaps
(
    const label index,
    const treeBoundBox& cubeBb
) const
{
    //return cubeBb.overlaps(allBb_[index]);

    //- Exact test of triangle intersecting bb

    // Quick rejection.
    if (!cubeBb.overlaps(allBb_[index]))
    {
        return false;
    }

    // Triangle points
    const pointField& points = surface_.points();
    const labelledTri& f = surface_[index];
    const point& p0 = points[f[0]];
    const point& p1 = points[f[1]];
    const point& p2 = points[f[2]];

    // Check if one or more triangle point inside
    if (cubeBb.contains(p0) || cubeBb.contains(p1) || cubeBb.contains(p2))
    {
        // One or more points inside
        return true;
    }

    // Now we have the difficult case: all points are outside but connecting
    // edges might go through cube. Use fast intersection of bounding box.

    return triangleFuncs::intersectBb(p0, p1, p2, cubeBb);
}


bool Foam::octreeDataTriSurface::contains
(
    const label,
    const point&
) const
{
    notImplemented
    (
        "octreeDataTriSurface::contains(const label, const point&)"
    );

    return false;
}


bool Foam::octreeDataTriSurface::intersects
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    if (mag(surface_.faceNormals()[index]) < VSMALL)
    {
        return false;
    }

    const pointField& points = surface_.points();

    const labelledTri& f = surface_[index];

    triPointRef tri(points[f[0]], points[f[1]], points[f[2]]);

    const vector dir(end - start);

    // Disable picking up intersections behind us.
    scalar oldTol = intersection::setPlanarTol(0.0);

    pointHit inter = tri.ray(start, dir, intersection::HALF_RAY);

    intersection::setPlanarTol(oldTol);

    if (inter.hit() && inter.distance() <= mag(dir))
    {
        // Note: no extra test on whether intersection is in front of us
        // since using half_ray AND zero tolerance. (note that tolerance
        // is used to look behind us)
        intersectionPoint = inter.hitPoint();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::octreeDataTriSurface::findTightest
(
    const label index,
    const point& sample,
    treeBoundBox& tightest
) const
{

    // get nearest and furthest away vertex
    point myNear, myFar;
    allBb_[index].calcExtremities(sample, myNear, myFar);

    const point dist = myFar - sample;
    scalar myFarDist = mag(dist);

    point tightestNear, tightestFar;
    tightest.calcExtremities(sample, tightestNear, tightestFar);

    scalar tightestFarDist = mag(tightestFar - sample);

    if (tightestFarDist < myFarDist)
    {
        // Keep current tightest.
        return false;
    }
    else
    {
        // Construct bb around sample and myFar
        const point dist2(fabs(dist.x()), fabs(dist.y()), fabs(dist.z()));

        tightest.min() = sample - dist2;
        tightest.max() = sample + dist2;

        return true;
    }
}


// Determine numerical value of sign of sample compared to shape at index
Foam::scalar Foam::octreeDataTriSurface::calcSign
(
    const label index,
    const point& sample,
    vector& n
) const
{
    n = surface_.faceNormals()[index];

    const labelledTri& tri = surface_[index];

    // take vector from sample to any point on triangle (we use vertex 0)
    vector vec = sample - surface_.points()[tri[0]];

    vec /= mag(vec) + VSMALL;

    return n & vec;
}


// Calculate nearest point to sample on/in shapei. !Does not set nearest
Foam::scalar Foam::octreeDataTriSurface::calcNearest
(
    const label index,
    const point& sample,
    point&
) const
{
    return mag(nearestPoint(index, sample) - sample);
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataTriSurface::calcNearest
(
    const label index,
    const linePointRef& ln,
    point& linePt,
    point& shapePt
) const
{
    notImplemented
    (
        "octreeDataTriSurface::calcNearest"
        "(const label, const linePointRef&, point& linePt, point&)"
    );
    return GREAT;
}


void Foam::octreeDataTriSurface::write
(
    Ostream& os,
    const label index
) const
{
    os << surface_[index] << token::SPACE << allBb_[index];
}


// ************************************************************************* //
