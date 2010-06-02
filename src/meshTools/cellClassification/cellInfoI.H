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

#include "cellClassification.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 information
inline bool Foam::cellInfo::update
(
    const cellInfo& w2,
    const label thisFaceI,
    const label thisCellI,
    const label neighbourFaceI,
    const label neighbourCellI
)
{
    if
    (
        (w2.type() == cellClassification::NOTSET)
     || (w2.type() == cellClassification::CUT)
    )
    {
        FatalErrorIn("cellInfo::update(const cellInfo&)")
            << "Problem: trying to propagate NOTSET or CUT type:" << w2.type()
            << " into cell/face with type:" << type() << endl
            << "thisFaceI:" << thisFaceI
            << "  thisCellI:" << thisCellI
            << "  neighbourFaceI:" << neighbourFaceI
            << "  neighbourCellI:" << neighbourCellI
            << abort(FatalError);
        return false;
    }

    if (type() == cellClassification::NOTSET)
    {
        type_ = w2.type();

        return true;
    }

    if (type() == cellClassification::CUT)
    {
        // Reached boundary. Stop.
        return false;
    }

    if (type() == w2.type())
    {
        // Should never happen; already checked in meshWave
        return false;
    }

    // Two conflicting types
    FatalErrorIn("cellInfo::update(const cellInfo&)")
        << "Problem: trying to propagate conflicting types:" << w2.type()
        << " into cell/face with type:" << type() << endl
        << "thisFaceI:" << thisFaceI
        << "  thisCellI:" << thisCellI
        << "  neighbourFaceI:" << neighbourFaceI
        << "  neighbourCellI:" << neighbourCellI
        << abort(FatalError);

    return false;
}
    

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::cellInfo::cellInfo()
:
    type_(cellClassification::NOTSET)
{}


// Construct from components
inline Foam::cellInfo::cellInfo(const label type)
:
    type_(type)
{}


// Construct as copy
inline Foam::cellInfo::cellInfo(const cellInfo& w2)
:
    type_(w2.type())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::cellInfo::valid() const
{
    return type_ != cellClassification::NOTSET;
}


// No geometric data so never any problem on cyclics
inline bool Foam::cellInfo::sameGeometry
(
    const polyMesh&,
    const cellInfo& w2,
    const scalar tol
)
 const
{
    return true;
}


// No geometric data.
inline void Foam::cellInfo::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// No geometric data.
inline void Foam::cellInfo::transform
(
    const polyMesh&,
    const tensor& rotTensor
)
{}


// No geometric data.
inline void Foam::cellInfo::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// Update this with neighbour information
inline bool Foam::cellInfo::updateCell
(
    const polyMesh&,
    const label thisCellI,
    const label neighbourFaceI,
    const cellInfo& neighbourInfo,
    const scalar tol
)
{
    return update
    (
        neighbourInfo,
        -1,
        thisCellI,
        neighbourFaceI,
        -1
    );
}    


// Update this with neighbour information
inline bool Foam::cellInfo::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const label neighbourCellI,
    const cellInfo& neighbourInfo,
    const scalar tol
)
{
    return update
    (
        neighbourInfo,
        thisFaceI,
        -1,
        -1,
        neighbourCellI
    );
}    

// Update this with neighbour information
inline bool Foam::cellInfo::updateFace
(
    const polyMesh&,
    const label thisFaceI,
    const cellInfo& neighbourInfo,
    const scalar tol
)
{
    return update
    (
        neighbourInfo,
        thisFaceI,
        -1,
        -1,
        -1
    );
}    


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::cellInfo::operator==(const Foam::cellInfo& rhs) const
{
    return type() == rhs.type();
}


inline bool Foam::cellInfo::operator!=(const Foam::cellInfo& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
