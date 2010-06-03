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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::sumData::sumData()
:
    oldFace_(-1),
    sum_(0.0),
    count_(0)
{}


// Construct from components
inline Foam::sumData::sumData
(
    const label oldFace,
    const scalar sum,
    const label count
)
:
    oldFace_(oldFace),
    sum_(sum),
    count_(count)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::sumData::valid() const
{
    return oldFace_ != -1;
}


// No geometric data so never any problem on cyclics
inline bool Foam::sumData::sameGeometry
(
    const polyMesh&,
    const sumData&,
    const scalar
) const
{
    return true;
}


// No geometric data.
inline void Foam::sumData::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// No geometric data.
inline void Foam::sumData::transform
(
    const polyMesh&,
    const tensor& rotTensor
)
{}


// No geometric data.
inline void Foam::sumData::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{
    oldFace_ = -1;
}


// Update cell with neighbouring face information
inline bool Foam::sumData::updateCell
(
    const polyMesh&,
    const label thisCellI,
    const label neighbourFaceI,
    const sumData& neighbourInfo,
    const scalar tol
)
{
    if (!valid())
    {
        FatalErrorIn("sumData::updateCell") << "problem"
            << abort(FatalError);
        return false;
    }


    if (count_ == 0)
    {
        sum_ += neighbourInfo.sum();
        count_ = neighbourInfo.count() + 1;
        oldFace_ = neighbourFaceI;
        return true;
    }
    else
    {
        return false;
    }
}


// Update face with neighbouring cell information
inline bool Foam::sumData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const sumData& neighbourInfo,
    const scalar tol
)
{
    // From cell to its faces.

    // Check that face is opposite the previous one.
    const cell& cFaces = mesh.cells()[neighbourCellI];

    label wantedFaceI = cFaces.opposingFaceLabel
    (
        neighbourInfo.oldFace(),
        mesh.faces()
    );

    if (thisFaceI == wantedFaceI)
    {
        if (count_ != 0)
        {
            FatalErrorIn("sumData::updateFace") << "problem"
                << abort(FatalError);
            return false;
        }

        sum_ += neighbourInfo.sum();
        count_ = neighbourInfo.count();
        oldFace_ = thisFaceI;
        return true;
    }
    else
    {
        return false;
    }
}


// Update face with coupled face information
inline bool Foam::sumData::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const sumData& neighbourInfo,
    const scalar tol
)
{
    // From face to face (e.g. coupled faces)
    if (count_ == 0)
    {
        sum_ += neighbourInfo.sum();
        count_ = neighbourInfo.count();
        oldFace_ = thisFaceI;
        return true;
    }
    else
    {
        return false;
    }
}    


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::sumData::operator==(const Foam::sumData& rhs)
 const
{
    return
        oldFace() == rhs.oldFace()
     && sum() == rhs.sum()
     && count() == rhs.count();
}


inline bool Foam::sumData::operator!=(const Foam::sumData& rhs)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
