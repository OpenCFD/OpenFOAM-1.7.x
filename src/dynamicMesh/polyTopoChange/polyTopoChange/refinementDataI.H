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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::refinementData::refinementData()
:
    refinementCount_(-1),
    count_(-1)
{}


// Construct from components
inline Foam::refinementData::refinementData
(
    const label refinementCount,
    const label count
)
:
    refinementCount_(refinementCount),
    count_(count)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::refinementData::valid() const
{
    return count_ != -1;
}


// No geometric data so never any problem on cyclics
inline bool Foam::refinementData::sameGeometry
(
    const polyMesh&,
    const refinementData&,
    const scalar
) const
{
    return true;
}


// No geometric data.
inline void Foam::refinementData::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// No geometric data.
inline void Foam::refinementData::transform
(
    const polyMesh&,
    const tensor& rotTensor
)
{}


// No geometric data.
inline void Foam::refinementData::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// Update cell with neighbouring face information
inline bool Foam::refinementData::updateCell
(
    const polyMesh&,
    const label thisCellI,
    const label neighbourFaceI,
    const refinementData& neighbourInfo,
    const scalar tol
)
{
    if (!valid())
    {
        FatalErrorIn("refinementData::updateCell") << "problem"
            << abort(FatalError);
        return false;
    }


    // Check if more than 2:1 ratio. This is when I am not refined but neighbour
    // is and neighbour already had higher cell level.
    if
    (
        neighbourInfo.isRefined()
    && !isRefined()
    &&  neighbourInfo.refinementCount() > refinementCount()
    )
    {
        count_ = refinementCount();
        return true;
    }



    // Count from neighbour face by the time it reaches the current cell.
    label transportedFaceCount;

    if (neighbourInfo.isRefined())
    {
        // refined so passes through two cells.
        transportedFaceCount = max(0, neighbourInfo.count()-2);
    }
    else
    {
        // unrefined.
        transportedFaceCount = max(0, neighbourInfo.count()-1);
    }

    if (count_ >= transportedFaceCount)
    {
        return false;
    }
    else
    {
        count_ = transportedFaceCount;

        return true;
    }
}


// Update face with neighbouring cell information
inline bool Foam::refinementData::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const label neighbourCellI,
    const refinementData& neighbourInfo,
    const scalar tol
)
{
    // From cell to its faces.
    if (!valid())
    {
        refinementCount_ = neighbourInfo.refinementCount();
        count_ = neighbourInfo.count();

        return true;
    }

    if (count_ >= neighbourInfo.count())
    {
        return false;
    }
    else
    {
        refinementCount_ = neighbourInfo.refinementCount();
        count_ = neighbourInfo.count();

        return true;
    }
}


// Update face with coupled face information
inline bool Foam::refinementData::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const refinementData& neighbourInfo,
    const scalar tol
)
{
    // From face to face (e.g. coupled faces)
    if (!valid())
    {
        refinementCount_ = neighbourInfo.refinementCount();
        count_ = neighbourInfo.count();

        return true;
    }

    if (count_ >= neighbourInfo.count())
    {
        return false;
    }
    else
    {
        refinementCount_ = neighbourInfo.refinementCount();
        count_ = neighbourInfo.count();

        return true;
    }
}    


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::refinementData::operator==(const Foam::refinementData& rhs)
 const
{
    return count() == rhs.count() && refinementCount() == rhs.refinementCount();
}


inline bool Foam::refinementData::operator!=(const Foam::refinementData& rhs)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
