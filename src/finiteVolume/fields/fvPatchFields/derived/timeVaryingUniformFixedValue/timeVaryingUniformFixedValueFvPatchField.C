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

#include "timeVaryingUniformFixedValueFvPatchField.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeSeries_()
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeSeries_(dict)
{
   if (dict.found("value"))
   {
       fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
   }
   else
   {
       updateCoeffs();
   }
}


template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_)
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    timeSeries_(ptf.timeSeries_)
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    timeSeries_(ptf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    fvPatchField<Type>::operator==
    (
        timeSeries_(this->db().time().timeOutputValue())
    );
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingUniformFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    timeSeries_.write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
