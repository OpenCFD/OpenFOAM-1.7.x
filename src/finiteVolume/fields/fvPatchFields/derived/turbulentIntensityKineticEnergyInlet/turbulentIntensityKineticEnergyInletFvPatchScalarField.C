/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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

#include "turbulentIntensityKineticEnergyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    intensity_(0.05)
{}

Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const turbulentIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    intensity_(ptf.intensity_)
{}

Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    intensity_(readScalar(dict.lookup("intensity")))
{
    if (intensity_ < 0 || intensity_ > 1)
    {
        FatalErrorIn
        (
            "turbulentIntensityKineticEnergyInletFvPatchScalarField::"
            "turbulentIntensityKineticEnergyInletFvPatchScalarField"
            "(const fvPatch& p, const DimensionedField<scalar, volMesh>& iF, "
            "const dictionary& dict)"
        )   << "Turbulence intensity should be specified as a fraction 0-1 "
               "of the mean velocity\n"
               "    value given is " << intensity_
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}

Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const turbulentIntensityKineticEnergyInletFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    intensity_(ptf.intensity_)
{}


Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
    const turbulentIntensityKineticEnergyInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    intensity_(ptf.intensity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    operator==(1.5*sqr(intensity_)*magSqr(Up));

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("intensity") << intensity_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        turbulentIntensityKineticEnergyInletFvPatchScalarField
    );
}

// ************************************************************************* //
