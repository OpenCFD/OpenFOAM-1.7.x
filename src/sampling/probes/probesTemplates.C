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

#include "probes.H"
#include "volFields.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- comparison operator for probes class
template<class T>
class isNotEqOp
{
public:

    void operator()(T& x, const T& y) const
    {
        const T unsetVal(-VGREAT*pTraits<T>::one);

        if (x != unsetVal)
        {
            // Keep x.

            // Note:chould check for y != unsetVal but multiple sample cells
            // already handled in read().
        }
        else
        {
            // x is not set. y might be.
            x = y;
        }
    }
};

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::label Foam::probes::countFields
(
    fieldGroup<Type>& fieldList,
    const wordList& fieldTypes
) const
{
    fieldList.setSize(fieldNames_.size());
    label nFields = 0;

    forAll(fieldNames_, fieldI)
    {
        if
        (
            fieldTypes[fieldI]
         == GeometricField<Type, fvPatchField, volMesh>::typeName
        )
        {
            fieldList[nFields] = fieldNames_[fieldI];
            nFields++;
        }
    }

    fieldList.setSize(nFields);

    return nFields;
}


template<class Type>
void Foam::probes::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    Field<Type> values = sample(vField);

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& probeStream = *probeFilePtrs_[vField.name()];

        probeStream << setw(w) << vField.time().value();

        forAll(values, probeI)
        {
            probeStream << ' ' << setw(w) << values[probeI];
        }
        probeStream << endl;
    }
}


template <class Type>
void Foam::probes::sampleAndWrite
(
    const fieldGroup<Type>& fields
)
{
    forAll(fields, fieldI)
    {
        if (loadFromFiles_)
        {
            sampleAndWrite
            (
                GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        fields[fieldI],
                        obr_.time().timeName(),
                        refCast<const polyMesh>(obr_),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    refCast<const fvMesh>(obr_)
                )
            );
        }
        else
        {
            objectRegistry::const_iterator iter = obr_.find(fields[fieldI]);

            if
            (
                iter != obr_.end()
             && iter()->type()
             == GeometricField<Type, fvPatchField, volMesh>::typeName
            )
            {
                sampleAndWrite
                (
                    obr_.lookupObject
                    <GeometricField<Type, fvPatchField, volMesh> >
                    (
                        fields[fieldI]
                    )
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::probes::sample
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    tmp<Field<Type> > tValues
    (
        new Field<Type>(probeLocations_.size(), unsetVal)
    );

    Field<Type>& values = tValues();

    forAll(probeLocations_, probeI)
    {
        if (cellList_[probeI] >= 0)
        {
            values[probeI] = vField[cellList_[probeI]];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::probes::sample(const word& fieldName) const
{
    return sample
    (
        obr_.lookupObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            fieldName
        )
    );
}


// ************************************************************************* //
