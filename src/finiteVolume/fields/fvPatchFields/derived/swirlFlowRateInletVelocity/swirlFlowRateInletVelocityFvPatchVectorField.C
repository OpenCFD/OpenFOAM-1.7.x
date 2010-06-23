/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "swirlFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(0),
    phiName_("phi"),
    rhoName_("rho"),
    rpm_(0)
{}


Foam::
swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_)
{}


Foam::
swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    flowRate_(readScalar(dict.lookup("flowRate"))),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rpm_(readScalar(dict.lookup("rpm")))
{}


Foam::
swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_)
{}


Foam::
swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swirlFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar totArea   = gSum(patch().magSf());
    // a simpler way of doing this would be nice
    const scalar avgU = -flowRate_/totArea;

    const vector avgCenter = gSum(patch().Cf()*patch().magSf())/totArea;
    const vector avgNormal = gSum(patch().Sf())/totArea;

    // Update angular velocity - convert [rpm] to [rad/s]
    vectorField tangentialVelocity =
    (
        (rpm_*mathematicalConstant::pi/30.0)
      * (patch().Cf() - avgCenter) ^ avgNormal
    );

    tmp<vectorField> n = patch().nf();

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // volumetric flow-rate
        operator==(tangentialVelocity + n*avgU);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // mass flow-rate
        operator==(tangentialVelocity + n*avgU/rhop);
    }
    else
    {
        FatalErrorIn
        (
            "swirlFlowRateInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << nl << exit(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::swirlFlowRateInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("flowRate") << flowRate_ << token::END_STATEMENT << nl;
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    os.writeKeyword("rpm") << rpm_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       swirlFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
