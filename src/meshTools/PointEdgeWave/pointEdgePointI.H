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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
inline bool Foam::pointEdgePoint::update
(
    const point& pt,
    const pointEdgePoint& w2,
    const scalar tol
)
{
    scalar dist2 = magSqr(pt - w2.origin());

    if (!valid())
    {
        // current not yet set so use any value
        distSqr_ = dist2;
        origin_ = w2.origin();

        return true;
    }        

    scalar diff = distSqr_ - dist2;

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < SMALL) || ((distSqr_ > SMALL) && (diff/distSqr_ < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr_ = dist2;
        origin_ = w2.origin();

        return true;
    }
}

    
// Update this with w2 (information on same point)
inline bool Foam::pointEdgePoint::update
(
    const pointEdgePoint& w2,
    const scalar tol
)
{
    if (!valid())
    {
        // current not yet set so use any value
        distSqr_ = w2.distSqr();
        origin_ = w2.origin();

        return true;
    }        

    scalar diff = distSqr_ - w2.distSqr();

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < SMALL) || ((distSqr_ > SMALL) && (diff/distSqr_ < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr_ =  w2.distSqr();
        origin_ = w2.origin();

        return true;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::pointEdgePoint::pointEdgePoint()
:
    origin_(greatPoint),
    distSqr_(GREAT)
{}


// Construct from origin, distance
inline Foam::pointEdgePoint::pointEdgePoint
(
    const point& origin, 
    const scalar distSqr
)
:
    origin_(origin),
    distSqr_(distSqr)
{}


// Construct as copy
inline Foam::pointEdgePoint::pointEdgePoint(const pointEdgePoint& wpt)
:
    origin_(wpt.origin()),
    distSqr_(wpt.distSqr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::pointEdgePoint::origin() const
{
    return origin_;
}


inline Foam::scalar Foam::pointEdgePoint::distSqr() const
{
    return distSqr_;
}


inline bool Foam::pointEdgePoint::valid() const
{
    return origin_ != greatPoint;
}


// Checks for cyclic points
inline bool Foam::pointEdgePoint::sameGeometry
(
    const pointEdgePoint& w2,
    const scalar tol
) const
{
    scalar diff = Foam::mag(distSqr() - w2.distSqr());

    if (diff < SMALL)
    {
        return true;
    }
    else
    {
        if ((distSqr() > SMALL) && ((diff/distSqr()) < tol))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


inline void Foam::pointEdgePoint::leaveDomain
(
    const polyPatch& patch,
    const label patchPointI,
    const point& coord
)
{
    origin_ -= coord;
}


inline void Foam::pointEdgePoint::transform(const tensor& rotTensor)
{
    origin_ = Foam::transform(rotTensor, origin_);
}


// Update absolute geometric quantities. Note that distance (distSqr_)
// is not affected by leaving/entering domain.
inline void Foam::pointEdgePoint::enterDomain
(
    const polyPatch& patch,
    const label patchPointI,
    const point& coord
)
{
    // back to absolute form
    origin_ += coord;
}


// Update this with information from connected edge
inline bool Foam::pointEdgePoint::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const label edgeI,
    const pointEdgePoint& edgeInfo,
    const scalar tol
)
{
    return 
        update
        (
            mesh.points()[pointI],
            edgeInfo,
            tol
        );
    }    


// Update this with new information on same point
inline bool Foam::pointEdgePoint::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const pointEdgePoint& newPointInfo,
    const scalar tol
)
{
    return
        update
        (
            mesh.points()[pointI],
            newPointInfo,
            tol
        );
}    


// Update this with new information on same point. No extra information.
inline bool Foam::pointEdgePoint::updatePoint
(
    const pointEdgePoint& newPointInfo,
    const scalar tol
)
{
    return update(newPointInfo, tol);
}    


// Update this with information from connected point
inline bool Foam::pointEdgePoint::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointI,
    const pointEdgePoint& pointInfo,
    const scalar tol
)
{
    const pointField& points = mesh.points();

    const edge& e = mesh.edges()[edgeI];

    const point edgeMid(0.5*(points[e[0]] + points[e[1]]));

    return
        update
        (
            edgeMid,
            pointInfo,
            tol
        );
}    


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointEdgePoint::operator==(const Foam::pointEdgePoint& rhs)
 const
{
    return origin() == rhs.origin();
}


inline bool Foam::pointEdgePoint::operator!=(const Foam::pointEdgePoint& rhs)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
