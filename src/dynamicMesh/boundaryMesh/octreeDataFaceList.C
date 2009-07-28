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

\*---------------------------------------------------------------------------*/

#include "octreeDataFaceList.H"
#include "octree.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::octreeDataFaceList, 0);

Foam::scalar Foam::octreeDataFaceList::tol = 1E-6;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::octreeDataFaceList::calcBb()
{
    allBb_.setSize(faceLabels_.size());
    allBb_ = treeBoundBox
    (
        vector(GREAT, GREAT, GREAT),
        vector(-GREAT, -GREAT, -GREAT)
    );

    forAll (faceLabels_, faceLabelI)
    {
        label faceI = faceLabels_[faceLabelI];

        // Update bb of face
        treeBoundBox& myBb = allBb_[faceLabelI];

        const face& f = mesh_.localFaces()[faceI];

        forAll(f, fp)
        {
            const point& coord = mesh_.localPoints()[f[fp]];

            myBb.min() = min(myBb.min(), coord);
            myBb.max() = max(myBb.max(), coord);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from all faces in bMesh
Foam::octreeDataFaceList::octreeDataFaceList(const bMesh& mesh)
:
    mesh_(mesh),
    faceLabels_(mesh_.size()),
    allBb_(mesh_.size())
{
    forAll(faceLabels_, faceI)
    {
        faceLabels_[faceI] = faceI;
    }

    // Generate tight fitting bounding box
    calcBb();
}


// Construct from selected faces in bMesh
Foam::octreeDataFaceList::octreeDataFaceList
(
    const bMesh& mesh,
    const labelList& faceLabels
)
:
    mesh_(mesh),
    faceLabels_(faceLabels),
    allBb_(faceLabels.size())
{
    // Generate tight fitting bounding box
    calcBb();
}



// Construct as copy
Foam::octreeDataFaceList::octreeDataFaceList(const octreeDataFaceList& shapes)
:
    mesh_(shapes.mesh()),
    faceLabels_(shapes.faceLabels()),
    allBb_(shapes.allBb())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::octreeDataFaceList::~octreeDataFaceList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::label Foam::octreeDataFaceList::getSampleType
(
    const octree<octreeDataFaceList>& oc,
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
    treeBoundBox tightest(treeBoundBox::greatBox);

    scalar tightestDist = GREAT;

    label index = oc.findNearest(sample, tightest, tightestDist);

    if (index == -1)
    {
        FatalErrorIn
        (
            "octreeDataFaceList::getSampleType"
            "(octree<octreeDataFaceList>&, const point&)"
        )   << "Could not find " << sample << " in octree."
            << abort(FatalError);
    }

    label faceI = faceLabels_[index];

    // Get actual intersection point on face

    if (debug & 2)
    {
        Pout<< "getSampleType : sample:" << sample
            << " nearest face:" << faceI;
    }

    const face& f = mesh_.localFaces()[faceI];

    const pointField& points = mesh_.localPoints();

    pointHit curHit = f.nearestPoint(sample, points);

    //
    // 1] Check whether sample is above face
    //

    if (curHit.hit())
    {
        // Simple case. Compare to face normal.

        if (debug & 2)
        {
            Pout<< " -> face hit:" << curHit.hitPoint()
                << " comparing to face normal " << mesh_.faceNormals()[faceI]
                << endl;
        }
        return octree<octreeDataFaceList>::getVolType
        (
            mesh_.faceNormals()[faceI],
            sample - curHit.hitPoint()
        );
    }

    if (debug & 2)
    {
        Pout<< " -> face miss:" << curHit.missPoint();
    }

    //
    // 2] Check whether intersection is on one of the face vertices or
    //    face centre
    //

    // typical dimension as sqrt of face area.
    scalar typDim = sqrt(mag(f.normal(points))) + VSMALL;

    forAll(f, fp)
    {
        if ((mag(points[f[fp]] - curHit.missPoint())/typDim) < tol)
        {
            // Face intersection point equals face vertex fp

            if (debug & 2)
            {
                    Pout<< " -> face point hit :" << points[f[fp]]
                        << " point normal:" << mesh_.pointNormals()[f[fp]]
                        << " distance:"
                        << mag(points[f[fp]] - curHit.missPoint())/typDim
                        << endl;
            }
            return octree<octreeDataFaceList>::getVolType
            (
                mesh_.pointNormals()[f[fp]],
                sample - curHit.missPoint()
            );
        }
    }

    // Get face centre
    point ctr(f.centre(points));

    if ((mag(ctr - curHit.missPoint())/typDim) < tol)
    {
        // Face intersection point equals face centre. Normal at face centre
        // is already average of face normals

        if (debug & 2)
        {
            Pout<< " -> centre hit:" << ctr
                << " distance:"
                << mag(ctr - curHit.missPoint())/typDim
                << endl;
        }

        return octree<octreeDataFaceList>::getVolType
        (
            mesh_.faceNormals()[faceI],
            sample - curHit.missPoint()
        );
    }


    //
    // 3] Get the 'real' edge the face intersection is on
    //

    const labelList& myEdges = mesh_.faceEdges()[faceI];

    forAll(myEdges, myEdgeI)
    {
        const edge& e = mesh_.edges()[myEdges[myEdgeI]];

        pointHit edgeHit =
            line<point, const point&>
            (
                points[e.start()],
                points[e.end()]
            ).nearestDist(sample);

        point edgePoint;
        if (edgeHit.hit())
        {
            edgePoint = edgeHit.hitPoint();
        }
        else
        {
            edgePoint = edgeHit.missPoint();
        }


        if ((mag(edgePoint - curHit.missPoint())/typDim) < tol)
        {
            // Face intersection point lies on edge e

            // Calculate edge normal (wrong: uses face normals instead of
            // triangle normals)
            const labelList& myFaces = mesh_.edgeFaces()[myEdges[myEdgeI]];

            vector edgeNormal(vector::zero);

            forAll(myFaces, myFaceI)
            {
                edgeNormal += mesh_.faceNormals()[myFaces[myFaceI]];
            }

            if (debug & 2)
            {
                Pout<< " -> real edge hit point:" << edgePoint
                    << " comparing to edge normal:" << edgeNormal
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return octree<octreeDataFaceList>::getVolType
            (
                edgeNormal,
                sample - curHit.missPoint()
            );
        }
    }


    //
    // 4] Get the internal edge (vertex - faceCentre) the face intersection
    //    is on
    //

    forAll(f, fp)
    {
        pointHit edgeHit =
            line<point, const point&>
            (
                points[f[fp]],
                ctr
            ).nearestDist(sample);

        point edgePoint;
        if (edgeHit.hit())
        {
            edgePoint = edgeHit.hitPoint();
        }
        else
        {
            edgePoint = edgeHit.missPoint();
        }

        if ((mag(edgePoint - curHit.missPoint())/typDim) < tol)
        {
            // Face intersection point lies on edge between two face triangles

            // Calculate edge normal as average of the two triangle normals
            label fpPrev = f.rcIndex(fp);
            label fpNext = f.fcIndex(fp);

            vector e = points[f[fp]] - ctr;
            vector ePrev = points[f[fpPrev]] - ctr;
            vector eNext = points[f[fpNext]] - ctr;

            vector nLeft = ePrev ^ e;
            nLeft /= mag(nLeft) + VSMALL;

            vector nRight = e ^ eNext;
            nRight /= mag(nRight) + VSMALL;

            if (debug & 2)
            {
                Pout<< " -> internal edge hit point:" << edgePoint
                    << " comparing to edge normal "
                    << 0.5*(nLeft + nRight)
                    << endl;
            }

            // Found face intersection point on this edge. Compare to edge
            // normal
            return octree<octreeDataFaceList>::getVolType
            (
                0.5*(nLeft + nRight),
                sample - curHit.missPoint()
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

    return octree<octreeDataFaceList>::UNKNOWN;
}


bool Foam::octreeDataFaceList::overlaps
(
    const label index,
    const treeBoundBox& sampleBb
) const
{
    return sampleBb.overlaps(allBb_[index]);
}


bool Foam::octreeDataFaceList::contains
(
    const label,
    const point&
) const
{
    notImplemented
    (
        "octreeDataFaceList::contains(const label, const point&)"
    );
    return false;
}


bool Foam::octreeDataFaceList::intersects
(
    const label index,
    const point& start,
    const point& end,
    point& intersectionPoint
) const
{
    label faceI = faceLabels_[index];

    const face& f = mesh_.localFaces()[faceI];

    const vector dir(end - start);

    // Disable picking up intersections behind us.
    scalar oldTol = intersection::setPlanarTol(0.0);

    pointHit inter =
        f.ray
        (
            start,
            dir,
            mesh_.localPoints(),
            intersection::HALF_RAY,
            intersection::VECTOR
        );

    intersection::setPlanarTol(oldTol);

    if (inter.hit() && inter.distance() <= mag(dir))
    {
        intersectionPoint = inter.hitPoint();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::octreeDataFaceList::findTightest
(
    const label index,
    const point& sample,
    treeBoundBox& tightest
) const
{
    // Get nearest and furthest away vertex
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
Foam::scalar Foam::octreeDataFaceList::calcSign
(
    const label index,
    const point& sample,
    vector&
) const
{
    label faceI = faceLabels_[index];

    const face& f = mesh_.localFaces()[faceI];

    point ctr = f.centre(mesh_.localPoints());

    vector vec = sample - ctr;

    vec /= mag(vec) + VSMALL;

    return mesh_.faceNormals()[faceI] & vec;
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataFaceList::calcNearest
(
    const label index,
    const point& sample,
    point& nearest
) const
{
    label faceI = faceLabels_[index];

    const face& f = mesh_.localFaces()[faceI];

    pointHit nearHit = f.nearestPoint(sample, mesh_.localPoints());

    if (nearHit.hit())
    {
        nearest = nearHit.hitPoint();
    }
    else
    {
        nearest = nearHit.missPoint();
    }

    if (debug & 1)
    {
        point ctr = f.centre(mesh_.localPoints());

        scalar sign = mesh_.faceNormals()[faceI] & (sample - nearest);

        Pout<< "octreeDataFaceList::calcNearest : "
            << "sample:" << sample
            << "  faceI:" << faceI
            << "  ctr:" << ctr
            << "  sign:" << sign
            << "  nearest point:" << nearest
            << "  distance to face:" << nearHit.distance()
            << endl;
    }
    return nearHit.distance();
}


void Foam::octreeDataFaceList::write
(
    Ostream& os,
    const label index
) const
{
    os << faceLabels_[index] << " " << allBb_[index];
}


// ************************************************************************* //
