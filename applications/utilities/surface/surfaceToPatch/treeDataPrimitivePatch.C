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

\*---------------------------------------------------------------------------*/

#include "treeDataPrimitivePatch.H"
#include "indexedOctree.H"
#include "triangleFuncs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::scalar
Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
tolSqr = sqr(1E-6);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::treeBoundBox
Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
calcBb
(
    const pointField& points,
    const face& f
)
{
    treeBoundBox bb(points[f[0]], points[f[0]]);

    for (label fp = 1; fp < f.size(); fp++)
    {
        const point& p = points[f[fp]];

        bb.min() = min(bb.min(), p);
        bb.max() = max(bb.max(), p);
    }
    return bb;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
update()
{
    if (cacheBb_)
    {
        bbs_.setSize(patch_.size());

        forAll(patch_, i)
        {
            bbs_[i] = calcBb(patch_.points(), patch_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
treeDataPrimitivePatch
(
    const bool cacheBb,
    const PrimitivePatch<Face, FaceList, PointField, PointType>& patch
)
:
    patch_(patch),
    cacheBb_(cacheBb)
{
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::pointField
Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
points() const
{
    pointField cc(patch_.size());

    forAll(patch_, i)
    {
        cc[i] = patch_[i].centre(patch_.points());
    }

    return cc;
}


//- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
//  Only makes sense for closed surfaces.
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::label
Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
getVolumeType
(
    const indexedOctree
    <
        treeDataPrimitivePatch
        <
            Face,
            FaceList,
            PointField,
            PointType
        >
    >& oc,
    const point& sample
) const
{
    // Need to determine whether sample is 'inside' or 'outside'
    // Done by finding nearest face. This gives back a face which is
    // guaranteed to contain nearest point. This point can be
    // - in interior of face: compare to face normal
    // - on edge of face: compare to edge normal
    // - on point of face: compare to point normal
    // Unfortunately the octree does not give us back the intersection point
    // or where on the face it has hit so we have to recreate all that
    // information.


    // Find nearest face to sample
    pointIndexHit info = oc.findNearest(sample, sqr(GREAT));

    if (info.index() == -1)
    {
        FatalErrorIn
        (
            "treeDataPrimitivePatch::getSampleType"
            "(indexedOctree<treeDataPrimitivePatch>&, const point&)"
        )   << "Could not find " << sample << " in octree."
            << abort(FatalError);
    }


    // Get actual intersection point on face
    label faceI = info.index();

    if (debug & 2)
    {
        Pout<< "getSampleType : sample:" << sample
            << " nearest face:" << faceI;
    }

    const pointField& points = patch_.localPoints();
    const face& f = patch_.localFaces()[faceI];

    // Retest to classify where on face info is. Note: could be improved. We
    // already have point.

    pointHit curHit = f.nearestPoint(sample, points);
    const vector area = f.normal(points);
    const point& curPt = curHit.rawPoint();

    //
    // 1] Check whether sample is above face
    //

    if (curHit.hit())
    {
        // Nearest point inside face. Compare to face normal.

        if (debug & 2)
        {
            Pout<< " -> face hit:" << curPt
                << " comparing to face normal " << area << endl;
        }
        return indexedOctree<treeDataPrimitivePatch>::getSide
        (
            area,
            sample - curPt
        );
    }

    if (debug & 2)
    {
        Pout<< " -> face miss:" << curPt;
    }

    //
    // 2] Check whether intersection is on one of the face vertices or
    //    face centre
    //

    const scalar typDimSqr = mag(area) + VSMALL;

    forAll(f, fp)
    {
        if ((magSqr(points[f[fp]] - curPt)/typDimSqr) < tolSqr)
        {
            // Face intersection point equals face vertex fp

            // Calculate point normal (wrong: uses face normals instead of
            // triangle normals)

            return indexedOctree<treeDataPrimitivePatch>::getSide
            (
                patch_.pointNormals()[f[fp]],
                sample - curPt
            );
        }
    }

    const point fc(f.centre(points));

    if ((magSqr(fc - curPt)/typDimSqr) < tolSqr)
    {
        // Face intersection point equals face centre. Normal at face centre
        // is already average of face normals

        if (debug & 2)
        {
            Pout<< " -> centre hit:" << fc
                << " distance:" << magSqr(fc - curPt)/typDimSqr << endl;
        }

        return indexedOctree<treeDataPrimitivePatch>::getSide
        (
            area,
            sample - curPt
        );
    }



    //
    // 3] Get the 'real' edge the face intersection is on
    //

    const labelList& fEdges = patch_.faceEdges()[faceI];

    forAll(fEdges, fEdgeI)
    {
        label edgeI = fEdges[fEdgeI];
        const edge& e = patch_.edges()[edgeI];

        pointHit edgeHit = e.line(points).nearestDist(sample);

        if ((magSqr(edgeHit.rawPoint() - curPt)/typDimSqr) < tolSqr)
        {
            // Face intersection point lies on edge e

            // Calculate edge normal (wrong: uses face normals instead of
            // triangle normals)
            const labelList& eFaces = patch_.edgeFaces()[edgeI];

            vector edgeNormal(vector::zero);

            forAll(eFaces, i)
            {
                edgeNormal += patch_.faceNormal()[eFaces[i]];
            }

            if (debug & 2)
            {
                Pout<< " -> real edge hit point:" << edgeHit.rawPoint()
                    << " comparing to edge normal:" << edgeNormal
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return indexedOctree<treeDataPrimitivePatch>::getSide
            (
                edgeNormal,
                sample - curPt
            );
        }
    }


    //
    // 4] Get the internal edge the face intersection is on
    //

    forAll(f, fp)
    {
        pointHit edgeHit = linePointRef
        (
            points[f[fp]],
            fc
        ).nearestDist(sample);

        if ((magSqr(edgeHit.rawPoint() - curPt)/typDimSqr) < tolSqr)
        {
            // Face intersection point lies on edge between two face triangles

            // Calculate edge normal as average of the two triangle normals
            vector e = points[f[fp]] - fc;
            vector ePrev = points[f[f.rcIndex(fp)]] - fc;
            vector eNext = points[f[f.fcIndex(fp)]] - fc;

            vector nLeft = ePrev ^ e;
            nLeft /= mag(nLeft) + VSMALL;

            vector nRight = e ^ eNext;
            nRight /= mag(nRight) + VSMALL;

            if (debug & 2)
            {
                Pout<< " -> internal edge hit point:" << edgeHit.rawPoint()
                    << " comparing to edge normal "
                    << 0.5*(nLeft + nRight)
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return indexedOctree<treeDataPrimitivePatch>::getSide
            (
                0.5*(nLeft + nRight),
                sample - curPt
            );
        }
    }

    if (debug & 2)
    {
        Pout<< "Did not find sample " << sample
            << " anywhere related to nearest face " << faceI << endl
            << "Face:";

        forAll(f, fp)
        {
            Pout<< "    vertex:" << f[fp] << "  coord:" << points[f[fp]]
                << endl;
        }
    }

    // Can't determine status of sample with respect to nearest face.
    // Either
    // - tolerances are wrong. (if e.g. face has zero area)
    // - or (more likely) surface is not closed.

    return indexedOctree<treeDataPrimitivePatch>::UNKNOWN;
}


// Check if any point on shape is inside cubeBb.
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
bool
Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
overlaps
(
    const label index,
    const treeBoundBox& cubeBb
) const
{
    // 1. Quick rejection: bb does not intersect face bb at all
    if (cacheBb_)
    {
        if (!cubeBb.overlaps(bbs_[index]))
        {
            return false;
        }
    }
    else
    {
        if (!cubeBb.overlaps(calcBb(patch_.points(), patch_[index])))
        {
            return false;
        }
    }


    // 2. Check if one or more face points inside

    const pointField& points = patch_.points();
    const face& f = patch_[index];

    forAll(f, fp)
    {
        if (cubeBb.contains(points[f[fp]]))
        {
            return true;
        }
    }

    // 3. Difficult case: all points are outside but connecting edges might
    // go through cube. Use triangle-bounding box intersection.
    const point fc = f.centre(points);

    forAll(f, fp)
    {
        bool triIntersects = triangleFuncs::intersectBb
        (
            points[f[fp]],
            points[f[f.fcIndex(fp)]],
            fc,
            cubeBb
        );

        if (triIntersects)
        {
            return true;
        }
    }
    return false;
}


// Calculate nearest point to sample. Updates (if any) nearestDistSqr, minIndex,
// nearestPoint.
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
findNearest
(
    const labelList& indices,
    const point& sample,

    scalar& nearestDistSqr,
    label& minIndex,
    point& nearestPoint
) const
{
    const pointField& points = patch_.points();

    forAll(indices, i)
    {
        label index = indices[i];

        const face& f = patch_[index];

        pointHit nearHit = f.nearestPoint(sample, points);
        scalar distSqr = sqr(nearHit.distance());

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            minIndex = index;
            nearestPoint = nearHit.rawPoint();
        }
    }
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
bool
Foam::treeDataPrimitivePatch<Face, FaceList, PointField, PointType>::
intersects
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    // Do quick rejection test
    if (cacheBb_)
    {
        const treeBoundBox& faceBb = bbs_[index];

        if ((faceBb.posBits(start) & faceBb.posBits(end)) != 0)
        {
            // start and end in same block outside of faceBb.
            return false;
        }
    }

    const pointField& points = patch_.points();
    const face& f = patch_[index];
    const point fc = f.centre(points);
    const vector dir(end - start);

    pointHit inter = patch_[index].intersection
    (
        start,
        dir,
        fc,
        points,
        intersection::HALF_RAY
    );

    if (inter.hit() && inter.distance() <= 1)
    {
        // Note: no extra test on whether intersection is in front of us
        // since using half_ray
        intersectionPoint = inter.hitPoint();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
