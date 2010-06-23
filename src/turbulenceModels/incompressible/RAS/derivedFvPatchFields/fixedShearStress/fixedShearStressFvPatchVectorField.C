/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fixedShearStressFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "RASModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    tau0_(vector::zero)
{}


fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    tau0_(dict.lookupOrDefault<vector>("tau", vector::zero))
{
    fvPatchField<vector>::operator=(patchInternalField());
}


fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    tau0_(ptf.tau0_)
{}


fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    tau0_(ptf.tau0_)
{}


fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    tau0_(ptf.tau0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedShearStressFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];

    const vectorField Ui = Uw.patchInternalField();

    vector tauHat = tau0_/(mag(tau0_) + ROOTVSMALL);

    const scalarField& ry = patch().deltaCoeffs();

    scalarField nuEffw = rasModel.nuEff()().boundaryField()[patchI];

    vectorField UwUpdated =
        tauHat*(tauHat & (tau0_*(1.0/(ry*nuEffw)) + Ui));

    operator==(UwUpdated);

    if (debug)
    {
        vectorField nHat = this->patch().nf();
        volSymmTensorField Reff = rasModel.devReff();
        Info << "tau : " << (nHat & Reff.boundaryField()[patchI])() << endl;
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void fixedShearStressFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    os.writeKeyword("tau") << tau0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedShearStressFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
