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

#include "directMappedVelocityFluxFixedValueFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "directMappedPatchBase.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("undefinedPhi")
{}


directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const directMappedVelocityFluxFixedValueFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "directMappedVelocityFluxFixedValueFvPatchField::"
            "directMappedVelocityFluxFixedValueFvPatchField\n"
            "(\n"
            "    const directMappedVelocityFluxFixedValueFvPatchField&,\n"
            "    const fvPatch&,\n"
            "    const DimensionedField<vector, volMesh>&,\n"
            "    const fvPatchFieldMapper&\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.lookup("phi"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "directMappedVelocityFluxFixedValueFvPatchField::"
            "directMappedVelocityFluxFixedValueFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const directMappedVelocityFluxFixedValueFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    phiName_(ptf.phiName_)
{}


directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const directMappedVelocityFluxFixedValueFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void directMappedVelocityFluxFixedValueFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the directMappedPatchBase
    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    (
        directMappedVelocityFluxFixedValueFvPatchField::patch().patch()
    );
    const mapDistribute& distMap = mpp.map();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const word& fieldName = dimensionedInternalField().name();
    const volVectorField& UField = nbrMesh.lookupObject<volVectorField>
    (
        fieldName
    );

    surfaceScalarField& phiField = const_cast<surfaceScalarField&>
    (
        nbrMesh.lookupObject<surfaceScalarField>(phiName_)
    );

    vectorField newUValues;
    scalarField newPhiValues;

    switch (mpp.mode())
    {
        case directMappedPolyPatch::NEARESTFACE:
        {
            vectorField allUValues(nbrMesh.nFaces(), vector::zero);
            scalarField allPhiValues(nbrMesh.nFaces(), 0.0);

            forAll(UField.boundaryField(), patchI)
            {
                const fvPatchVectorField& Upf = UField.boundaryField()[patchI];
                const scalarField& phipf = phiField.boundaryField()[patchI];

                label faceStart = Upf.patch().patch().start();

                forAll(Upf, faceI)
                {
                    allUValues[faceStart + faceI] = Upf[faceI];
                    allPhiValues[faceStart + faceI] = phipf[faceI];
                }
            }

            mapDistribute::distribute
            (
                Pstream::defaultCommsType,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                distMap.constructMap(),
                allUValues
            );
            newUValues = patch().patchSlice(allUValues);

            mapDistribute::distribute
            (
                Pstream::defaultCommsType,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                distMap.constructMap(),
                allPhiValues
            );
            newPhiValues = patch().patchSlice(allPhiValues);

            break;
        }
        case directMappedPolyPatch::NEARESTPATCHFACE:
        {
            const label nbrPatchID = nbrMesh.boundaryMesh().findPatchID
            (
                mpp.samplePatch()
            );

            newUValues = UField.boundaryField()[nbrPatchID];

            mapDistribute::distribute
            (
                Pstream::defaultCommsType,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                distMap.constructMap(),
                newUValues
            );

            newPhiValues = phiField.boundaryField()[nbrPatchID];

            mapDistribute::distribute
            (
                Pstream::defaultCommsType,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                distMap.constructMap(),
                newPhiValues
            );

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "directMappedVelocityFluxFixedValueFvPatchField::updateCoeffs()"
            )<< "patch can only be used in NEARESTPATCHFACE or NEARESTFACE "
             << "mode" << nl << abort(FatalError);
        }
    }

    operator==(newUValues);
    phiField.boundaryField()[patch().index()] == newPhiValues;

    fixedValueFvPatchVectorField::updateCoeffs();
}


void directMappedVelocityFluxFixedValueFvPatchField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    directMappedVelocityFluxFixedValueFvPatchField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
