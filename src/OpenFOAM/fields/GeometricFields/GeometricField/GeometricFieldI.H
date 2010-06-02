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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
inline const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::null()
{
    return *reinterpret_cast< GeometricField<Type, PatchField, GeoMesh>* >(0);
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline
const typename
Foam::GeometricField<Type, PatchField, GeoMesh>::DimensionedInternalField&
Foam::GeometricField<Type, PatchField, GeoMesh>::
dimensionedInternalField() const
{
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline
const typename
Foam::GeometricField<Type, PatchField, GeoMesh>::InternalField&
Foam::GeometricField<Type, PatchField, GeoMesh>::internalField() const
{
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline const typename Foam::GeometricField<Type, PatchField, GeoMesh>::
GeometricBoundaryField&
Foam::GeometricField<Type, PatchField, GeoMesh>::boundaryField() const
{
    return boundaryField_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline Foam::label
Foam::GeometricField<Type, PatchField, GeoMesh>::timeIndex() const
{
    return timeIndex_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline Foam::label&
Foam::GeometricField<Type, PatchField, GeoMesh>::timeIndex()
{
    return timeIndex_;
}


// ************************************************************************* //
