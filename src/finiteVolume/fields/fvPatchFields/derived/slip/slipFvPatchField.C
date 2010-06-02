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

#include "slipFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
slipFvPatchField<Type>::slipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
   basicSymmetryFvPatchField<Type>(p, iF)
{}


template<class Type>
slipFvPatchField<Type>::slipFvPatchField
(
    const slipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basicSymmetryFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
slipFvPatchField<Type>::slipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    basicSymmetryFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
slipFvPatchField<Type>::slipFvPatchField
(
    const slipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    basicSymmetryFvPatchField<Type>(ptf, iF)
{}


template<class Type>
slipFvPatchField<Type>::slipFvPatchField
(
    const slipFvPatchField<Type>& ptf
)
:
    basicSymmetryFvPatchField<Type>(ptf)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
