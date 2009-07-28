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

#include "directMappedFixedValueFvPatchField.H"
#include "directMappedPatchBase.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    setAverage_(false),
    average_(pTraits<Type>::zero)
{}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const directMappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "directMappedFixedValueFvPatchField<Type>::"
            "directMappedFixedValueFvPatchField\n"
            "(\n"
            "    const directMappedFixedValueFvPatchField<Type>&,\n"
            "    const fvPatch&,\n"
            "    const Field<Type>&,\n"
            "    const fvPatchFieldMapper&\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    setAverage_(readBool(dict.lookup("setAverage"))),
    average_(pTraits<Type>(dict.lookup("average")))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "directMappedFixedValueFvPatchField<Type>::"
            "directMappedFixedValueFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    // Force calculation of schedule (uses parallel comms)
    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    (
        this->patch().patch()
    );
    (void)mpp.map().schedule();
}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const directMappedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const directMappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void directMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Get the scheduling information from the directMappedPatchBase
    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    (
        directMappedFixedValueFvPatchField<Type>::patch().patch()
    );
    const mapDistribute& distMap = mpp.map();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const word& fldName = this->dimensionedInternalField().name();

    // Result of obtaining remote values
    Field<Type> newValues;

    switch (mpp.mode())
    {
        case directMappedPatchBase::NEARESTCELL:
        {
            if (mpp.sameRegion())
            {
                newValues = this->internalField();
            }
            else
            {
                newValues = nbrMesh.lookupObject<fieldType>
                (
                    fldName
                ).internalField();
            }
            mapDistribute::distribute
            (
                Pstream::defaultCommsType,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                distMap.constructMap(),
                newValues
            );

            break;
        }
        case directMappedPatchBase::NEARESTPATCHFACE:
        {
            const label nbrPatchID = nbrMesh.boundaryMesh().findPatchID
            (
                mpp.samplePatch()
            );
            if (nbrPatchID < 0)
            {
                FatalErrorIn
                (
                    "void directMappedFixedValueFvPatchField<Type>::"
                    "updateCoeffs()"
                )<< "Unable to find sample patch " << mpp.samplePatch()
                 << " in region " << mpp.sampleRegion()
                 << " for patch " << this->patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = nbrMesh.lookupObject<fieldType>
            (
                fldName
            );
            newValues = nbrField.boundaryField()[nbrPatchID];
            mapDistribute::distribute
            (
                Pstream::defaultCommsType,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                distMap.constructMap(),
                newValues
            );

            break;
        }
        case directMappedPatchBase::NEARESTFACE:
        {
            Field<Type> allValues(nbrMesh.nFaces(), pTraits<Type>::zero);

            const fieldType& nbrField = nbrMesh.lookupObject<fieldType>
            (
                fldName
            );
            forAll(nbrField.boundaryField(), patchI)
            {
                const fvPatchField<Type>& pf =
                    nbrField.boundaryField()[patchI];
                label faceStart = pf.patch().patch().start();

                forAll(pf, faceI)
                {
                    allValues[faceStart++] = pf[faceI];
                }
            }

            mapDistribute::distribute
            (
                Pstream::defaultCommsType,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                distMap.constructMap(),
                allValues
            );

            newValues = this->patch().patchSlice(allValues);

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "directMappedFixedValueFvPatchField<Type>::updateCoeffs()"
            )<< "Unknown sampling mode: " << mpp.mode()
             << nl << abort(FatalError);
        }
    }

    if (setAverage_)
    {
        Type averagePsi =
            gSum(this->patch().magSf()*newValues)
           /gSum(this->patch().magSf());

        if (mag(averagePsi)/mag(average_) > 0.5)
        {
            newValues *= mag(average_)/mag(averagePsi);
        }
        else
        {
            newValues += (average_ - averagePsi);
        }
    }

    this->operator==(newValues);

    if (debug)
    {
        Info<< "directMapped on field:" << fldName
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void directMappedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;
    os.writeKeyword("average") << average_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
