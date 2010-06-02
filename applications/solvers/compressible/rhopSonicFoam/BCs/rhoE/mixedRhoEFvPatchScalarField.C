/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "mixedRhoEFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixedRhoEFvPatchScalarField::mixedRhoEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


mixedRhoEFvPatchScalarField::mixedRhoEFvPatchScalarField
(
    const mixedRhoEFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


mixedRhoEFvPatchScalarField::mixedRhoEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


mixedRhoEFvPatchScalarField::mixedRhoEFvPatchScalarField
(
    const mixedRhoEFvPatchScalarField& ptpsf
)
:
    mixedFvPatchScalarField(ptpsf)
{}


mixedRhoEFvPatchScalarField::mixedRhoEFvPatchScalarField
(
    const mixedRhoEFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mixedRhoEFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void mixedRhoEFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


void mixedRhoEFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& rhop =
        patch().lookupPatchField<volScalarField, scalar>("rho");
    const fvPatchField<vector>& rhoUp =
        patch().lookupPatchField<volVectorField, vector>("rhoU");
//    fvPatchField<scalar>& Tp =
//        patch().lookupPatchField<volScalarField, scalar>("T");

    const volScalarField& T = db().lookupObject<volScalarField>("T");
    const label patchi = patch().index();
    fvPatchScalarField& Tp = 
        const_cast<fvPatchScalarField&>(T.boundaryField()[patchi]);

    Tp.evaluate();

    const dictionary& thermodynamicProperties = db().lookupObject<IOdictionary>
    (
        "thermodynamicProperties"
    );

    dimensionedScalar Cv(thermodynamicProperties.lookup("Cv"));

    valueFraction() = rhop.snGrad()/
        (rhop.snGrad() - rhop*this->patch().deltaCoeffs());

    refValue() = 0.5*rhop*magSqr(rhoUp/rhop);
    refGrad() =
        rhop*Cv.value()*Tp.snGrad()
      + (
            refValue() 
          - (0.5*rhop.patchInternalField()*
                magSqr(rhoUp.patchInternalField()/rhop.patchInternalField()))
        )*patch().deltaCoeffs();

    mixedFvPatchScalarField::updateCoeffs();
}


void mixedRhoEFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("valueFraction") << valueFraction()
        << token::END_STATEMENT << endl;
    os.writeKeyword("refValue") << refValue() << token::END_STATEMENT << endl;
    os.writeKeyword("refGrad") << refGrad() << token::END_STATEMENT << endl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, mixedRhoEFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
