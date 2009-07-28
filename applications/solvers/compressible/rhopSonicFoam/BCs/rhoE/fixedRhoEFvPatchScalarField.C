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

#include "fixedRhoEFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedRhoEFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& thermodynamicProperties = db().lookupObject<IOdictionary>
    (
        "thermodynamicProperties"
    );

    dimensionedScalar Cv(thermodynamicProperties.lookup("Cv"));

    const fvPatchScalarField& rhop =
        patch().lookupPatchField<volScalarField, scalar>("rho");

    const fvPatchVectorField& rhoUp =
        patch().lookupPatchField<volVectorField, vector>("rhoU");

    const fvPatchScalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>("T");

    operator==(rhop*(Cv.value()*Tp + 0.5*magSqr(rhoUp/rhop)));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void fixedRhoEFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedRhoEFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
