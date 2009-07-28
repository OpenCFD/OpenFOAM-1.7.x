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

#include "expDirectionMixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
expDirectionMixedFvPatchField<Type>::expDirectionMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{}


template<class Type>
expDirectionMixedFvPatchField<Type>::expDirectionMixedFvPatchField
(
    const expDirectionMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    refGrad_(ptf.refGrad_, mapper),
    valueFraction_(ptf.valueFraction_, mapper)
{}


template<class Type>
expDirectionMixedFvPatchField<Type>::expDirectionMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    refValue_("refValue", dict, p.size()),
    refGrad_("refGradient", dict, p.size()),
    valueFraction_("valueFraction", dict, p.size())
{
    evaluate();
}


template<class Type>
expDirectionMixedFvPatchField<Type>::expDirectionMixedFvPatchField
(
    const expDirectionMixedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


template<class Type>
expDirectionMixedFvPatchField<Type>::expDirectionMixedFvPatchField
(
    const expDirectionMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void expDirectionMixedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    refValue_.autoMap(m);
    refGrad_.autoMap(m);
    valueFraction_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void expDirectionMixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const expDirectionMixedFvPatchField<Type>& edmptf =
        refCast<const expDirectionMixedFvPatchField<Type> >(ptf);

    refValue_.rmap(edmptf.refValue_, addr);
    refGrad_.rmap(edmptf.refGrad_, addr);
    valueFraction_.rmap(edmptf.valueFraction_, addr);
}


template<class Type>
tmp<Field<Type> > expDirectionMixedFvPatchField<Type>::snGrad() const
{
    const vectorField& nHat = patch().faceNormals();

    Field<Type> gradValue =
        patchInternalField() + refGrad_/patch().deltaCoeffs();

    Field<Type> mixedValue =
        nHat*(nHat & refValue_)
      + gradValue - nHat*(nHat & gradValue);

    return
        valueFraction_*
            (mixedValue - patchInternalField())*patch().deltaCoeffs()
      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
void expDirectionMixedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated())
    {
        updateCoeffs();
    }

    const vectorField& nHat = patch().faceNormals();

    Field<Type> gradValue =
        patchInternalField() + refGrad_/patch().deltaCoeffs();

    Field<Type> mixedValue =
        nHat*(nHat & refValue_)
      + gradValue - nHat*(nHat & gradValue);

    Field<Type>::operator=
    (
        valueFraction_*mixedValue + (1.0 - valueFraction_)*gradValue
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > expDirectionMixedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one)*(1.0 - valueFraction_);
}


template<class Type>
tmp<Field<Type> > expDirectionMixedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    const vectorField& nHat = patch().faceNormals();

    Field<Type> gradValue =
        patchInternalField() + refGrad_/patch().deltaCoeffs();

    Field<Type> mixedValue =
        nHat*(nHat & refValue_)
      + gradValue - nHat*(nHat & gradValue);

    return
        valueFraction_*mixedValue
      + (1.0 - valueFraction_)*refGrad_/patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > expDirectionMixedFvPatchField<Type>::
gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*valueFraction_*patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > expDirectionMixedFvPatchField<Type>::
gradientBoundaryCoeffs() const
{
    const vectorField& nHat = patch().faceNormals();

    Field<Type> gradValue =
        patchInternalField() + refGrad_/patch().deltaCoeffs();

    Field<Type> mixedValue =
        nHat*(nHat & refValue_)
      + gradValue - nHat*(nHat & gradValue);

    return
        valueFraction_*patch().deltaCoeffs()*mixedValue
      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
void expDirectionMixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
