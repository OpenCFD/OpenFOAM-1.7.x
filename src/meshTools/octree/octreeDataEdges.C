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

#include "octreeDataEdges.H"

#include "line.H"
#include "labelList.H"
#include "octree.H"
#include "linePointRef.H"
#include "pointHit.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::octreeDataEdges, 0);

Foam::scalar Foam::octreeDataEdges::tol(1E-6);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from selected edges. Bounding box calculated.
Foam::octreeDataEdges::octreeDataEdges
(
    const edgeList& edges,
    const pointField& points,
    const labelList& edgeLabels
)
:
    edges_(edges),
    points_(points),
    edgeLabels_(edgeLabels),
    allBb_(edgeLabels_.size())
{
    // Generate tight fitting bounding box
    forAll(edgeLabels_, i)
    {
        label edgeI = edgeLabels_[i];

        const edge& e = edges_[edgeI];

        const point& a = points_[e.start()];
        const point& b = points_[e.end()];

        allBb_[i].min() = min(a, b);
        allBb_[i].max() = max(a, b);
    }
}


// Construct as copy
Foam::octreeDataEdges::octreeDataEdges(const octreeDataEdges& shapes)
:
    edges_(shapes.edges()),
    points_(shapes.points()),
    edgeLabels_(shapes.edgeLabels()),
    allBb_(shapes.allBb())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::octreeDataEdges::~octreeDataEdges()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::octreeDataEdges::getSampleType
(
    const octree<octreeDataEdges>&,
    const point&
) const
{
    return octree<octreeDataEdges>::UNKNOWN;
}


bool Foam::octreeDataEdges::overlaps
(
    const label index,
    const treeBoundBox& sampleBb
) const
{
    return sampleBb.overlaps(allBb_[index]);
}


bool Foam::octreeDataEdges::contains
(
    const label,
    const point&
) const
{
    notImplemented
    (
        "octreeDataEdges::contains(const label, const point&)"
    );
    return false;
}


bool Foam::octreeDataEdges::intersects
(
    const label,
    const point&,
    const point&,
    point&
) const
{
    notImplemented
    (
        "octreeDataEdges::intersects(const label, const point&"
        ", const point&, point&)"
    );
    return false;
}


bool Foam::octreeDataEdges::findTightest
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
Foam::scalar Foam::octreeDataEdges::calcSign
(
    const label,
    const point&,
    point& n
) const
{
    n = vector::zero;

    return 1;
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataEdges::calcNearest
(
    const label index,
    const point& sample,
    point& nearest
) const
{
    const edge& e = edges_[edgeLabels_[index]];

    pointHit nearHit = e.line(points_).nearestDist(sample);

    nearest = nearHit.rawPoint();

    return nearHit.distance();
}


// Calculate nearest point on/in shapei
Foam::scalar Foam::octreeDataEdges::calcNearest
(
    const label index,
    const linePointRef& sampleLine,
    point& sampleLinePt,
    point& shapePt
) const
{
    const edge& e = edges_[edgeLabels_[index]];

    linePointRef edgeLine(e.line(points_));

    return edgeLine.nearestDist(sampleLine, shapePt, sampleLinePt);
}
    

void Foam::octreeDataEdges::write
(
    Ostream& os,
    const label index
) const
{
    os << edgeLabels_[index] << " " << allBb_[index];
}


// ************************************************************************* //
