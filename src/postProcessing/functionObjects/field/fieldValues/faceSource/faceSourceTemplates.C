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

#include "faceSource.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "IOList.H"
#include "ListListOps.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::fieldValues::faceSource::setFieldValues
(
    const word& fieldName,
    List<Type>& values
) const
{
    values.setSize(faceId_.size(), pTraits<Type>::zero);

    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (obr_.foundObject<sf>(fieldName))
    {
        const sf& field = obr_.lookupObject<sf>(fieldName);

        forAll(values, i)
        {
            label faceI = faceId_[i];
            label patchI = facePatchId_[i];
            if (patchI >= 0)
            {
                values[i] = field.boundaryField()[patchI][faceI];
            }
            else
            {
                values[i] = field[faceI];
            }

            values[i] *= flipMap_[i];
        }

        return true;
    }
    else if (obr_.foundObject<vf>(fieldName))
    {
        const vf& field = obr_.lookupObject<vf>(fieldName);

        forAll(values, i)
        {
            label faceI = faceId_[i];
            label patchI = facePatchId_[i];
            if (patchI >= 0)
            {
                values[i] = field.boundaryField()[patchI][faceI];
            }
            else
            {
                FatalErrorIn
                (
                    "fieldValues::faceSource::setFieldValues"
                    "("
                        "const word&, "
                        "List<Type>&"
                    ") const"
                )   << type() << " " << name_ << ": "
                    << sourceTypeNames_[source_] << "(" << sourceName_ << "):"
                    << nl
                    << "    Unable to process internal faces for volume field "
                    << fieldName << nl << abort(FatalError);
            }

            values[i] *= flipMap_[i];
        }

        return true;
    }

    return false;
}


template<class Type>
Type Foam::fieldValues::faceSource::processValues
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
        case opAreaAverage:
        {
            tmp<scalarField> magSf = filterField(mesh().magSf());
            result = sum(values*magSf())/sum(magSf());
            break;
        }
        case opAreaIntegrate:
        {
            result = sum(values*filterField(mesh().magSf()));
            break;
        }
        case opWeightedAverage:
        {
            if (mesh().foundObject<volScalarField>(weightFieldName_))
            {
                tmp<scalarField> wField =
                    filterField
                    (
                        mesh().lookupObject<volScalarField>(weightFieldName_)
                    );
               result = sum(values*wField())/sum(wField());
            }
            else if (mesh().foundObject<surfaceScalarField>(weightFieldName_))
            {
                tmp<scalarField> wField =
                    filterField
                    (
                        mesh().lookupObject<surfaceScalarField>
                        (
                            weightFieldName_
                        )
                    );
               result = sum(values*wField())/sum(wField());
            }
            else
            {
                FatalErrorIn
                (
                    "fieldValues::faceSource::processValues"
                    "("
                        "List<Type>&"
                    ") const"
                )   << type() << " " << name_ << ": "
                    << sourceTypeNames_[source_] << "(" << sourceName_ << "):"
                    << nl
                    << "    Weight field " << weightFieldName_
                    << " must be either a " << volScalarField::typeName
                    << " or " << surfaceScalarField::typeName << nl
                    << abort(FatalError);
            }
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
bool Foam::fieldValues::faceSource::writeValues(const word& fieldName)
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
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::faceSource::filterField
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    tmp<Field<Type> > tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues();

    forAll(values, i)
    {
        label faceI = faceId_[i];
        label patchI = facePatchId_[i];
        if (patchI >= 0)
        {
            values[i] = field.boundaryField()[patchI][faceI];
        }
        else
        {
            FatalErrorIn
            (
                "fieldValues::faceSource::filterField"
                "("
                    "const GeometricField<Type, fvPatchField, volMesh>&"
                ") const"
            )   << type() << " " << name_ << ": "
                << sourceTypeNames_[source_] << "(" << sourceName_ << "):"
                << nl
                << "    Unable to process internal faces for volume field "
                << field.name() << nl << abort(FatalError);
        }

        values[i] *= flipMap_[i];
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::faceSource::filterField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
) const
{
    tmp<Field<Type> > tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues();

    forAll(values, i)
    {
        label faceI = faceId_[i];
        label patchI = facePatchId_[i];
        if (patchI >= 0)
        {
            values[i] = field.boundaryField()[patchI][faceI];
        }
        else
        {
            values[i] = field[faceI];
        }

        values[i] *= flipMap_[i];
    }

    return tvalues;
}


// ************************************************************************* //

