/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cyclicPointPatchField.H"
#include "Swap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cyclicPointPatchField<Type>::cyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicPointPatch>(p))
{}


template<class Type>
cyclicPointPatchField<Type>::cyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    cyclicPatch_(refCast<const cyclicPointPatch>(p))
{
    if (!isType<cyclicPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicPointPatchField<Type>::cyclicPointPatchField\n"
            "(\n"
            "    const pointPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not cyclic type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
cyclicPointPatchField<Type>::cyclicPointPatchField
(
    const cyclicPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicPointPatch>(p))
{
    if (!isType<cyclicPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicPointPatchField<Type>::cyclicPointPatchField\n"
            "(\n"
            "    const cyclicPointPatchField<Type>& ptf,\n"
            "    const pointPatch& p,\n"
            "    const DimensionedField<Type, pointMesh>& iF,\n"
            "    const pointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
cyclicPointPatchField<Type>::cyclicPointPatchField
(
    const cyclicPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void cyclicPointPatchField<Type>::swapAdd(Field<Type>& pField) const
{
    Field<Type> pf(this->patchInternalField(pField));

    const edgeList& pairs = cyclicPatch_.transformPairs();

    if (doTransform())
    {
        forAll(pairs, pairi)
        {
            Type tmp = pf[pairs[pairi][0]];
            pf[pairs[pairi][0]] = transform(forwardT()[0], pf[pairs[pairi][1]]);
            pf[pairs[pairi][1]] = transform(reverseT()[0], tmp);
        }
    }
    else
    {
        forAll(pairs, pairi)
        {
            Swap(pf[pairs[pairi][0]], pf[pairs[pairi][1]]);
        }
    }

    addToInternalField(pField, pf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
