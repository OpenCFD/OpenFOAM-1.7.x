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

#include "globalPointPatchField.H"
#include "PstreamCombineReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
globalPointPatchField<Type>::globalPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    globalPointPatch_(refCast<const globalPointPatch>(p))
{}


template<class Type>
globalPointPatchField<Type>::globalPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    globalPointPatch_(refCast<const globalPointPatch>(p))
{
    if (!isType<globalPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "globalPointPatchField<Type>::globalPointPatchField\n"
            "(\n"
            " const pointPatch& p,\n"
            " const Field<Type>& field,\n"
            " const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index()
            << " not processorPoint type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
globalPointPatchField<Type>::globalPointPatchField
(
    const globalPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    globalPointPatch_(refCast<const globalPointPatch>(ptf.patch()))
{
    if (!isType<globalPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "globalPointPatchField<Type>::globalPointPatchField\n"
            "(\n"
            " const globalPointPatchField<Type>& ptf,\n"
            " const pointPatch& p,\n"
            " const DimensionedField<Type, pointMesh>& iF,\n"
            " const pointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
globalPointPatchField<Type>::globalPointPatchField
(
    const globalPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    globalPointPatch_(refCast<const globalPointPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class Type>
globalPointPatchField<Type>::~globalPointPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class Type>
void globalPointPatchField<Type>::swapAdd(Field<Type>& pField) const
{
    // Create the global list and insert local values
    if (globalPointPatch_.globalPointSize() > 0)
    {
        Field<Type> lpf = patchInternalField(pField);
        const labelList& addr = globalPointPatch_.sharedPointAddr();

        Field<Type> gpf
        (
            globalPointPatch_.globalPointSize(),
            pTraits<Type>::zero
        );

        forAll(addr, i)
        {
            gpf[addr[i]] += lpf[i];
        }

        combineReduce(gpf, plusEqOp<Field<Type> >());

        // Extract local data
        forAll (addr, i)
        {
            lpf[i] = gpf[addr[i]];
        }

        setInInternalField(pField, lpf);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
