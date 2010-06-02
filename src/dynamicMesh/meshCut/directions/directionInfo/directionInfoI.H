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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::directionInfo::directionInfo()
:
    index_(-3),
    n_(vector::zero)
{}


// Construct from components
inline Foam::directionInfo::directionInfo
(
    const label index,
    const vector& n
)
:
    index_(index),
    n_(n)
{}


// Construct as copy
inline Foam::directionInfo::directionInfo(const directionInfo& w2)
:
    index_(w2.index()),
    n_(w2.n())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::directionInfo::valid() const
{
    return index_ != -3;
}


// No geometric data so never any problem on cyclics
inline bool Foam::directionInfo::sameGeometry
(
    const polyMesh&,
    const directionInfo& w2,
    const scalar tol
)
 const
{
    return true;
}


// index_ is already offset in face
inline void Foam::directionInfo::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// index_ is offset in face on other side. So reverse it here.
// (Note: f[0] on other domain is connected to f[0] in this domain,
//        f[1]      ,,                         f[size-1]   ,, 
// etc.)
inline void Foam::directionInfo::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{
    if (index_ >= 0)
    {
        const face& f = patch[patchFaceI];

        index_ = (f.size() - index_) % f.size();
    }
}


// No geometric data.
inline void Foam::directionInfo::transform
(
    const polyMesh&,
    const tensor& rotTensor
)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::directionInfo::operator==(const Foam::directionInfo& rhs)
 const
{
    return  index() == rhs.index() && n() == rhs.n();
}


inline bool Foam::directionInfo::operator!=(const Foam::directionInfo& rhs)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
