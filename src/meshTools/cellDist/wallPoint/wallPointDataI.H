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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Update this with w2 if w2 nearer to pt.
template <class Type>
inline bool wallPointData<Type>::update
(
    const point& pt,
    const wallPointData<Type>& w2,
    const scalar tol
)
{
    scalar dist2 = magSqr(pt - w2.origin());

    if (!valid())
    {
        // current not yet set so use any value
        distSqr() = dist2;
        origin() = w2.origin();
        data_ = w2.data();

        return true;
    }        

    scalar diff = distSqr() - dist2;

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < SMALL) || ((distSqr() > SMALL) && (diff/distSqr() < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr() = dist2;
        origin() = w2.origin();
        data_ = w2.data();

        return true;
    }
}
    

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
template <class Type>
inline wallPointData<Type>::wallPointData()
:
    wallPoint(),
    data_()
{}


// Construct from components
template <class Type>
inline wallPointData<Type>::wallPointData
(
    const point& origin,
    const Type& data,
    const scalar distSqr
)
:
    wallPoint(origin, distSqr),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Type>
inline const Type& wallPointData<Type>::data() const
{
    return data_;
}


template <class Type>
inline Type& wallPointData<Type>::data()
{
    return data_;
}


// Update this with w2 if w2 nearer to pt.
template <class Type>
inline bool wallPointData<Type>::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label,
    const wallPointData<Type>& neighbourWallInfo,
    const scalar tol
)
{
    const vectorField& cellCentres = mesh.primitiveMesh::cellCentres();

    return update
    (
        cellCentres[thisCellI],
        neighbourWallInfo,
        tol
    );
}    


// Update this with w2 if w2 nearer to pt.
template <class Type>
inline bool wallPointData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label,
    const wallPointData<Type>& neighbourWallInfo,
    const scalar tol
)
{
    const vectorField& faceCentres = mesh.faceCentres();

    return update
    (
        faceCentres[thisFaceI],
        neighbourWallInfo,
        tol
    );
}    


// Update this with w2 if w2 nearer to pt.
template <class Type>
inline bool wallPointData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const wallPointData<Type>& neighbourWallInfo,
    const scalar tol
)
{
    const vectorField& faceCentres = mesh.faceCentres();

    return update
    (
        faceCentres[thisFaceI],
        neighbourWallInfo,
        tol
    );
}    


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
