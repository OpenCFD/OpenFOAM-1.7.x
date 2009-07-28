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

#include "transform.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Returns the wanted level
inline Foam::label Foam::refinementDistanceData::wantedLevel(const point& pt)
 const
{
    const scalar distSqr = magSqr(pt-origin_);

    // Get the size at the origin level
    scalar levelSize = level0Size_/(1<<originLevel_);

    scalar r = 0;

    for (label level = originLevel_; level >= 0; --level)
    {
        // Current range
        r += levelSize;

        // Check if our distance is within influence sphere
        if (sqr(r) > distSqr)
        {
            return level;
        }

        // Lower level will have double the size
        levelSize *= 2;
    }
    return 0;
}


inline bool Foam::refinementDistanceData::update
(
    const point& pos,
    const refinementDistanceData& neighbourInfo,
    const scalar tol
)
{
    if (!valid())
    {
        if (!neighbourInfo.valid())
        {
            FatalErrorIn("refinementDistanceData::update(..)")
                << "problem" << abort(FatalError);
        }
        operator=(neighbourInfo);
        return true;
    }

    // Determine wanted level at current position.
    label cellLevel = wantedLevel(pos);

    // Determine wanted level coming through the neighbour
    label nbrLevel = neighbourInfo.wantedLevel(pos);

    if (nbrLevel > cellLevel)
    {
        operator=(neighbourInfo);
        return true;
    }
    else if (nbrLevel == cellLevel)
    {
        scalar myDistSqr = magSqr(pos-origin_);
        scalar nbrDistSqr = magSqr(pos - neighbourInfo.origin());
        scalar diff = myDistSqr - nbrDistSqr;

        if (diff < 0)
        {
            // already nearest
            return false;
        }

        if ((diff < SMALL) || ((myDistSqr > SMALL) && (diff/myDistSqr < tol)))
        {
            // don't propagate small changes
            return false;
        }
        else
        {
            // update with new values
            operator=(neighbourInfo);
            return true;
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::refinementDistanceData::refinementDistanceData()
:
    level0Size_(-1)
{}


// Construct from components
inline Foam::refinementDistanceData::refinementDistanceData
(
    const scalar level0Size,
    const point& origin,
    const label originLevel
)
:
    level0Size_(level0Size),
    origin_(origin),
    originLevel_(originLevel)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::refinementDistanceData::valid() const
{
    return level0Size_ != -1;
}


// No geometric data so never any problem on cyclics
inline bool Foam::refinementDistanceData::sameGeometry
(
    const polyMesh&,
    const refinementDistanceData&,
    const scalar
) const
{
    return true;
}


inline void Foam::refinementDistanceData::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{
    origin_ -= faceCentre;
}


inline void Foam::refinementDistanceData::transform
(
    const polyMesh&,
    const tensor& rotTensor
)
{
    origin_ = Foam::transform(rotTensor, origin_);
}


// Update absolute geometric quantities.
inline void Foam::refinementDistanceData::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{
    // back to absolute form
    origin_ += faceCentre;
}


// Update cell with neighbouring face information
inline bool Foam::refinementDistanceData::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const refinementDistanceData& neighbourInfo,
    const scalar tol
)
{
    const point& pos = mesh.cellCentres()[thisCellI];

    return update(pos, neighbourInfo, tol);
}


// Update face with neighbouring cell information
inline bool Foam::refinementDistanceData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const refinementDistanceData& neighbourInfo,
    const scalar tol
)
{
    const point& pos = mesh.faceCentres()[thisFaceI];

    return update(pos, neighbourInfo, tol);
}


// Update face with coupled face information
inline bool Foam::refinementDistanceData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const refinementDistanceData& neighbourInfo,
    const scalar tol
)
{
    const point& pos = mesh.faceCentres()[thisFaceI];

    return update(pos, neighbourInfo, tol);
}    


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::refinementDistanceData::operator==
(
    const Foam::refinementDistanceData& rhs
)
 const
{
    if (!valid())
    {
        if (!rhs.valid())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return
            level0Size_ == rhs.level0Size_
         && origin_ == rhs.origin_
         && originLevel_ == rhs.originLevel_;
    }
}


inline bool Foam::refinementDistanceData::operator!=
(
    const Foam::refinementDistanceData& rhs
)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
