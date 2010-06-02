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

#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if not yet set.
inline bool Foam::wallNormalInfo::update(const wallNormalInfo& w2)
{
    if (!w2.valid())
    {
        FatalErrorIn
        (
            "wallNormalInfo::update(const wallNormalInfo&)"
        ) << "Problem: w2 is not valid" << abort(FatalError);

        return false;
    }
    else if (valid())
    {
        // our already set. Stop any transfer
        return false;
    }
    else
    {
        normal_ = w2.normal();

        return true;
    }
}
    

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::wallNormalInfo::wallNormalInfo()
:
    normal_(greatVector)
{}


// Construct from normal
inline Foam::wallNormalInfo::wallNormalInfo(const vector& normal)
:
    normal_(normal)
{}


// Construct as copy
inline Foam::wallNormalInfo::wallNormalInfo(const wallNormalInfo& wpt)
:
    normal_(wpt.normal())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::wallNormalInfo::normal() const
{
    return normal_;
}


inline bool Foam::wallNormalInfo::valid() const
{
    return normal_ != greatVector;
}


// No geometric data so never any problem on cyclics
inline bool Foam::wallNormalInfo::sameGeometry
(
    const polyMesh&,
    const wallNormalInfo& w2,
    const scalar tol
) const
{
    return true;
}


// No geometric data.
inline void Foam::wallNormalInfo::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// No geometric data.
inline void Foam::wallNormalInfo::transform
(
    const polyMesh&,
    const tensor& rotTensor
)
{}


// No geometric data.
inline void Foam::wallNormalInfo::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// Update this with w2 if w2 nearer to pt.
inline bool Foam::wallNormalInfo::updateCell
(
    const polyMesh&,
    const label thisCellI,
    const label neighbourFaceI,
    const wallNormalInfo& neighbourWallInfo,
    const scalar tol
)
{
    return update(neighbourWallInfo);
}


// Update this with w2 if w2 nearer to pt.
inline bool Foam::wallNormalInfo::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const label neighbourCellI,
    const wallNormalInfo& neighbourWallInfo,
    const scalar tol
)
{
    return update(neighbourWallInfo);
}


// Update this with w2 if w2 nearer to pt.
inline bool Foam::wallNormalInfo::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const wallNormalInfo& neighbourWallInfo,
    const scalar tol
)
{
    return update(neighbourWallInfo);
}    


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wallNormalInfo::operator==(const Foam::wallNormalInfo& rhs)
 const
{
    return normal() == rhs.normal();
}


inline bool Foam::wallNormalInfo::operator!=(const Foam::wallNormalInfo& rhs)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
