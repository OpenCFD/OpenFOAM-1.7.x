/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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

#include "turbulentMixingLengthDissipationRateInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    mixingLength_(0.001)
{}

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    mixingLength_(ptf.mixingLength_)
{}

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    mixingLength_(readScalar(dict.lookup("mixingLength")))
{}

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    mixingLength_(ptf.mixingLength_)
{}

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    mixingLength_(ptf.mixingLength_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentMixingLengthDissipationRateInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup Cmu corresponding to the turbulence model selected
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const scalar Cmu =
        rasModel.coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);

    const scalar Cmu75 = pow(Cmu, 0.75);

    const fvPatchField<scalar>& kp =
        patch().lookupPatchField<volScalarField, scalar>("k");

    operator==(Cmu75*kp*sqrt(kp)/mixingLength_);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void turbulentMixingLengthDissipationRateInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("mixingLength")
        << mixingLength_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentMixingLengthDissipationRateInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
