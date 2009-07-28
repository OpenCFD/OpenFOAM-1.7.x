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

InClass
    Foam

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    transformList.C

\*---------------------------------------------------------------------------*/

#ifndef transformList_H
#define transformList_H

#include "transform.H"
#include "List.H"
#include "Map.H"
#include "EdgeMap.H"
#include "tensorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Extend transform to work on list.
template<class T>
List<T> transform
(
    const tensor& rotTensor,
    const UList<T>& field
);

//- Apply transformation to list. Either single transformation tensor
//  or one tensor per element.
template<class T>
void transformList(const tensorField&, UList<T>&);

template<class T>
void transformList(const tensorField&, Map<T>&);

template<class T>
void transformList(const tensorField&, EdgeMap<T>&);


template<>
inline void transformList(const tensorField&, UList<label>&)
{}
template<>
inline void transformList(const tensorField&, Map<label>&)
{}
template<>
inline void transformList(const tensorField&, EdgeMap<label>&)
{}


template<>
inline void transformList(const tensorField&, UList<scalar>&)
{}
template<>
inline void transformList(const tensorField&, Map<scalar>&)
{}
template<>
inline void transformList(const tensorField&, EdgeMap<scalar>&)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "transformList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
