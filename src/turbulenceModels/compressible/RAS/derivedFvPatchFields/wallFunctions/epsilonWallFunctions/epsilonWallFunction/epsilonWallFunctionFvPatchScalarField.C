/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "epsilonWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void epsilonWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("epsilonWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    UName_("U"),
    kName_("k"),
    GName_("RASModel::G"),
    rhoName_("rho"),
    muName_("mu"),
    mutName_("mut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    GName_(ptf.GName_),
    rhoName_(ptf.rhoName_),
    muName_(ptf.muName_),
    mutName_(ptf.mutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    GName_(dict.lookupOrDefault<word>("G", "RASModel::G")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    muName_(dict.lookupOrDefault<word>("mu", "mu")),
    mutName_(dict.lookupOrDefault<word>("mut", "mut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(ewfpsf),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    rhoName_(ewfpsf.rhoName_),
    muName_(ewfpsf.muName_),
    mutName_(ewfpsf.mutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_)
{
    checkType();
}


epsilonWallFunctionFvPatchScalarField::epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(ewfpsf, iF),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    rhoName_(ewfpsf.rhoName_),
    muName_(ewfpsf.muName_),
    mutName_(ewfpsf.mutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalar Cmu75 = pow(Cmu_, 0.75);
    const scalar yPlusLam = rasModel.yPlusLam(kappa_, E_);

    const scalarField& y = rasModel.y()[patch().index()];

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));

    volScalarField& epsilon = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(dimensionedInternalField().name()));

    const volScalarField& k = db().lookupObject<volScalarField>(kName_);

    const scalarField& rhow =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    const scalarField& muw =
        patch().lookupPatchField<volScalarField, scalar>(muName_);

    const scalarField& mutw =
        patch().lookupPatchField<volScalarField, scalar>(mutName_);

    const fvPatchVectorField& Uw =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const scalarField magGradUw = mag(Uw.snGrad());

    // Set epsilon and G
    forAll(mutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        epsilon[faceCellI] = Cmu75*pow(k[faceCellI], 1.5)/(kappa_*y[faceI]);

        G[faceCellI] =
            (mutw[faceI] + muw[faceI])
           *magGradUw[faceI]
           *Cmu25*sqrt(k[faceCellI])
           /(kappa_*y[faceI]);
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    fixedInternalValueFvPatchField<scalar>::updateCoeffs();
}


void epsilonWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedInternalValueFvPatchField<scalar>::evaluate(commsType);
}


void epsilonWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "mu", "mu", muName_);
    writeEntryIfDifferent<word>(os, "mut", "mut", mutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
