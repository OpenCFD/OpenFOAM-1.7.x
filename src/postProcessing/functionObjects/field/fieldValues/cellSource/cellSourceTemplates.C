/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "cellSource.H"
#include "volFields.H"
#include "IOList.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::fieldValues::cellSource::validField(const word& fieldName) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (obr_.foundObject<vf>(fieldName))
    {
        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::cellSource::setFieldValues
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (obr_.foundObject<vf>(fieldName))
    {
        return filterField(obr_.lookupObject<vf>(fieldName));
    }

    return tmp<Field<Type> >(new Field<Type>(0.0));
}


template<class Type>
Type Foam::fieldValues::cellSource::processValues
(
    const Field<Type>& values,
    const scalarField& V,
    const scalarField& weightField
) const
{
    Type result = pTraits<Type>::zero;
    switch (operation_)
    {
        case opSum:
        {
            result = sum(values);
            break;
        }
        case opVolAverage:
        {
            result = sum(values*V)/sum(V);
            break;
        }
        case opVolIntegrate:
        {
            result = sum(values*V);
            break;
        }
        case opWeightedAverage:
        {
            result = sum(values*weightField)/sum(weightField);
            break;
        }
        default:
        {
            // Do nothing
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fieldValues::cellSource::writeValues(const word& fieldName)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values = combineFields(setFieldValues<Type>(fieldName));

        scalarField V = combineFields(filterField(mesh().V()));

        scalarField weightField =
            combineFields(setFieldValues<scalar>(weightFieldName_));

        if (Pstream::master())
        {
            Type result = processValues(values, V, weightField);

            if (valueOutput_)
            {
                IOList<Type>
                (
                    IOobject
                    (
                        fieldName + "_" + sourceTypeNames_[source_] + "-"
                            + sourceName_,
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    values
                ).write();
            }


            outputFilePtr_()<< tab << result;

            if (log_)
            {
                Info<< "    " << operationTypeNames_[operation_]
                    << "(" << sourceName_ << ") for " << fieldName
                    <<  " = " << result << endl;
            }
        }
    }

    return ok;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::cellSource::filterField
(
    const Field<Type>& field
) const
{
    return tmp<Field<Type> >(new Field<Type>(field, cellId_));
}


// ************************************************************************* //

