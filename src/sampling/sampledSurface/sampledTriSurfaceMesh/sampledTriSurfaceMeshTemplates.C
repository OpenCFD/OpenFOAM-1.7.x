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

#include "sampledTriSurfaceMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledTriSurfaceMesh::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // One value per face
    tmp<Field<Type> > tvalues(new Field<Type>(cellLabels_.size()));
    Field<Type>& values = tvalues();

    forAll(cellLabels_, triI)
    {
        values[triI] = vField[cellLabels_[triI]];
    }

    return tvalues;
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledTriSurfaceMesh::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    tmp<Field<Type> > tvalues(new Field<Type>(pointToFace_.size()));
    Field<Type>& values = tvalues();

    forAll(pointToFace_, pointI)
    {
        label triI = pointToFace_[pointI];
        label cellI = cellLabels_[triI];

        values[pointI] = interpolator.interpolate(points()[pointI], cellI);
    }

    return tvalues;
}


// ************************************************************************* //
