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

#include "nuSgsWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nuSgsWallFunctionFvPatchScalarField::
nuSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    nuName_("nu"),
    kappa_(0.41),
    E_(9.8)
{}


nuSgsWallFunctionFvPatchScalarField::
nuSgsWallFunctionFvPatchScalarField
(
    const nuSgsWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    nuName_(ptf.nuName_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


nuSgsWallFunctionFvPatchScalarField::
nuSgsWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{}


nuSgsWallFunctionFvPatchScalarField::
nuSgsWallFunctionFvPatchScalarField
(
    const nuSgsWallFunctionFvPatchScalarField& nwfpsf
)
:
    fixedValueFvPatchScalarField(nwfpsf),
    UName_(nwfpsf.UName_),
    nuName_(nwfpsf.nuName_),
    kappa_(nwfpsf.kappa_),
    E_(nwfpsf.E_)
{}


nuSgsWallFunctionFvPatchScalarField::
nuSgsWallFunctionFvPatchScalarField
(
    const nuSgsWallFunctionFvPatchScalarField& nwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(nwfpsf, iF),
    UName_(nwfpsf.UName_),
    nuName_(nwfpsf.nuName_),
    kappa_(nwfpsf.kappa_),
    E_(nwfpsf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nuSgsWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const scalarField& ry = patch().deltaCoeffs();

    const fvPatchVectorField& U =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    scalarField magUp = mag(U.patchInternalField() - U);

    const scalarField& nuw =
        patch().lookupPatchField<volScalarField, scalar>(nuName_);

    scalarField& nuSgsw = *this;

    scalarField magFaceGradU = mag(U.snGrad());

    forAll(nuSgsw, facei)
    {
        scalar magUpara = magUp[facei];

        scalar utau = sqrt((nuSgsw[facei] + nuw[facei])*magFaceGradU[facei]);

        if (utau > VSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUpara/utau, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - utau/(ry[facei]*nuw[facei])
                    + magUpara/utau
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    - 1.0/(ry[facei]*nuw[facei])
                    - magUpara/sqr(utau)
                    - 1/E_*kUu*fkUu/utau;

                scalar utauNew = utau - f/df;
                err = mag((utau - utauNew)/utau);
                utau = utauNew;

            } while (utau > VSMALL && err > 0.01 && ++iter < 10);

            nuSgsw[facei] =
                max(sqr(max(utau, 0))/magFaceGradU[facei] - nuw[facei], 0.0);
        }
        else
        {
            nuSgsw[facei] = 0;
        }
    }

    fixedValueFvPatchScalarField::evaluate();
}


void nuSgsWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nuSgsWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
