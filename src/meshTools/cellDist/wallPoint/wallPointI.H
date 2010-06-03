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
inline bool Foam::wallPoint::update
(
    const point& pt,
    const wallPoint& w2,
    const scalar tol
)
{
    //Already done in calling algorithm
    //if (w2.origin() == origin_)
    //{
    //    // Shortcut. Same input so same distance.
    //    return false;
    //}

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
    

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::wallPoint::wallPoint()
:
    origin_(greatPoint),
    distSqr_(GREAT)
{}


// Construct from origin, distance
inline Foam::wallPoint::wallPoint(const point& origin, const scalar distSqr)
:
    origin_(origin), distSqr_(distSqr)
{}


// Construct as copy
inline Foam::wallPoint::wallPoint(const wallPoint& wpt)
:
    origin_(wpt.origin()), distSqr_(wpt.distSqr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::wallPoint::origin() const
{
    return origin_;
}


inline Foam::point& Foam::wallPoint::origin()
{
    return origin_;
}


inline Foam::scalar Foam::wallPoint::distSqr() const
{
    return distSqr_;
}


inline Foam::scalar& Foam::wallPoint::distSqr()
{
    return distSqr_;
}


inline bool Foam::wallPoint::valid() const
{
    return origin_ != greatPoint;
}


// Checks for cyclic faces
inline bool Foam::wallPoint::sameGeometry
(
    const polyMesh&,
    const wallPoint& w2,
    const scalar tol
)
 const
{
    scalar diff = mag(distSqr() - w2.distSqr());

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


inline void Foam::wallPoint::leaveDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point& faceCentre
)
{
    origin_ -= faceCentre;
}


inline void Foam::wallPoint::transform
(
    const polyMesh&,
    const tensor& rotTensor
)
{
    origin_ = Foam::transform(rotTensor, origin_);
}


// Update absolute geometric quantities. Note that distance (distSqr_)
// is not affected by leaving/entering domain.
inline void Foam::wallPoint::enterDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point& faceCentre
)
{
    // back to absolute form
    origin_ += faceCentre;
}


// Update this with w2 if w2 nearer to pt.
inline bool Foam::wallPoint::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const wallPoint& neighbourWallInfo,
    const scalar tol
)
{
    return 
        update
        (
            mesh.cellCentres()[thisCellI],
            neighbourWallInfo,
            tol
        );
    }    


// Update this with w2 if w2 nearer to pt.
inline bool Foam::wallPoint::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const wallPoint& neighbourWallInfo,
    const scalar tol
)
{
    return
        update
        (
            mesh.faceCentres()[thisFaceI],
            neighbourWallInfo,
            tol
        );
}    

// Update this with w2 if w2 nearer to pt.
inline bool Foam::wallPoint::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const wallPoint& neighbourWallInfo,
    const scalar tol
)
{
    return
        update
        (
            mesh.faceCentres()[thisFaceI],
            neighbourWallInfo,
            tol
        );
}    


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wallPoint::operator==(const Foam::wallPoint& rhs) const
{
    return origin() == rhs.origin();
}


inline bool Foam::wallPoint::operator!=(const Foam::wallPoint& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
