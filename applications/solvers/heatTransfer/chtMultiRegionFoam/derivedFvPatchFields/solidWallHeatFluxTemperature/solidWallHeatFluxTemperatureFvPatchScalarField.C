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

#include "solidWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    KName_("undefined-K")
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const solidWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    q_(ptf.q_, mapper),
    KName_(ptf.KName_)
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    q_("q", dict, p.size()),
    KName_(dict.lookup("K"))
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const solidWallHeatFluxTemperatureFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    q_(tppsf.q_),
    KName_(tppsf.KName_)
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const solidWallHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    q_(tppsf.q_),
    KName_(tppsf.KName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    q_.autoMap(m);
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const solidWallHeatFluxTemperatureFvPatchScalarField& hfptf =
        refCast<const solidWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    q_.rmap(hfptf.q_, addr);
}


Foam::tmp<Foam::scalarField>
Foam::solidWallHeatFluxTemperatureFvPatchScalarField::K() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (mesh.objectRegistry::foundObject<volScalarField>(KName_))
    {
        return patch().lookupPatchField<volScalarField, scalar>(KName_);
    }
    else if (mesh.objectRegistry::foundObject<volSymmTensorField>(KName_))
    {
        const symmTensorField& KWall =
            patch().lookupPatchField<volSymmTensorField, scalar>(KName_);

        vectorField n = patch().nf();

        return n & KWall & n;
    }
    else
    {
        FatalErrorIn
        (
            "solidWallHeatFluxTemperatureFvPatchScalarField::K()"
            " const"
        )   << "Did not find field " << KName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'K' to a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    gradient() = q_/K();

    fixedGradientFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(K()*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heatFlux:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    q_.writeEntry("q", os);
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        solidWallHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
