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

#include "processorPointPatchField.H"
#include "transformField.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
processorPointPatchField<Type>::processorPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    procPatch_(refCast<const processorPointPatch>(p))
{}


template<class Type>
processorPointPatchField<Type>::processorPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorPointPatch>(p))
{}


template<class Type>
processorPointPatchField<Type>::processorPointPatchField
(
    const processorPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorPointPatch>(ptf.patch()))
{}


template<class Type>
processorPointPatchField<Type>::processorPointPatchField
(
    const processorPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorPointPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
processorPointPatchField<Type>::~processorPointPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void processorPointPatchField<Type>::initSwapAdd(Field<Type>& pField) const
{
    if (Pstream::parRun())
    {
        // Get internal field into my point order
        Field<Type> pf(this->patchInternalField(pField));

        OPstream::write
        (
            Pstream::blocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(pf.begin()),
            pf.byteSize()
        );
    }
}


template<class Type>
void processorPointPatchField<Type>::swapAdd(Field<Type>& pField) const
{
    if (Pstream::parRun())
    {
        Field<Type> pnf(this->size());

        IPstream::read
        (
            Pstream::blocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(pnf.begin()),
            pnf.byteSize()
        );

        if (doTransform())
        {
            const processorPolyPatch& ppp = procPatch_.procPolyPatch();
            const tensorField& forwardT = ppp.forwardT();

            if (forwardT.size() == 1)
            {
                transform(pnf, forwardT[0], pnf);
            }
            else
            {
                const labelList& nonGlobalPatchPoints =
                    procPatch_.nonGlobalPatchPoints();
                const labelListList& pointFaces = ppp.pointFaces();

                forAll(nonGlobalPatchPoints, pfi)
                {
                    pnf[pfi] = transform
                    (
                        forwardT[pointFaces[nonGlobalPatchPoints[pfi]][0]],
                        pnf[pfi]
                    );
                }
            }
        }

        addToInternalField(pField, pnf);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
