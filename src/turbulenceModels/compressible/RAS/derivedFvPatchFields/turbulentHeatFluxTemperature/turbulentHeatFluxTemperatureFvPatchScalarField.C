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

#include "turbulentHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    rhoName_("rho")
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    q_(ptf.q_, mapper),
    rhoName_(ptf.rhoName_)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_("q", dict, p.size()),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& thftpsf
)
:
    fixedGradientFvPatchScalarField(thftpsf),
    q_(thftpsf.q_),
    rhoName_(thftpsf.rhoName_)
{}


turbulentHeatFluxTemperatureFvPatchScalarField::
turbulentHeatFluxTemperatureFvPatchScalarField
(
    const turbulentHeatFluxTemperatureFvPatchScalarField& thftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(thftpsf, iF),
    q_(thftpsf.q_),
    rhoName_(thftpsf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    q_.autoMap(m);
}


void turbulentHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const turbulentHeatFluxTemperatureFvPatchScalarField& thftptf =
        refCast<const turbulentHeatFluxTemperatureFvPatchScalarField>
        (
            ptf
        );

    q_.rmap(thftptf.q_, addr);
}


void turbulentHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const scalarField alphaEffp = rasModel.alphaEff()().boundaryField()[patchI];

    const basicThermo& thermo =
        db().lookupObject<basicThermo>("thermophysicalProperties");

//    const scalarField& Tp = thermo.T().boundaryField()[patchI];
    const scalarField& Tp = *this;

    const scalarField Cpp = thermo.Cp(Tp, patchI);

    const scalarField& rhop =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const scalar Ap = gSum(patch().magSf());

    gradient() = q_/(Ap*rhop*Cpp*alphaEffp);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void turbulentHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    q_.writeEntry("q", os);
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    gradient().writeEntry("gradient", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentHeatFluxTemperatureFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //

