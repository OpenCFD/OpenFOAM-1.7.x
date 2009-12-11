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
#include "ListListOps.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::fieldValues::cellSource::setFieldValues
(
    const word& fieldName,
    List<Type>& values
) const
{
    values.setSize(cellId_.size(), pTraits<Type>::zero);

    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (obr_.foundObject<vf>(fieldName))
    {
        const vf& field = obr_.lookupObject<vf>(fieldName);

        values = UIndirectList<Type>(field, cellId_);

        return true;
    }

    return false;
}


template<class Type>
Type Foam::fieldValues::cellSource::processValues
(
    const List<Type>& values
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
            tmp<scalarField> V = filterField(mesh().V());
            result = sum(values*V())/sum(V());
            break;
        }
        case opVolIntegrate:
        {
            result = sum(values*filterField(mesh().V()));
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
    List<List<Type> > allValues(Pstream::nProcs());

    bool validField =
        setFieldValues<Type>(fieldName, allValues[Pstream::myProcNo()]);

    if (validField)
    {
        Pstream::gatherList(allValues);

        if (Pstream::master())
        {
            List<Type> values =
                ListListOps::combine<List<Type> >
                (
                    allValues,
                    accessOp<List<Type> >()
                );

            Type result = processValues(values);

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

    return validField;
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

