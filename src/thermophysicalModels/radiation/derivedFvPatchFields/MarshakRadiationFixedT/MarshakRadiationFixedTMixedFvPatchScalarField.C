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

#include "MarshakRadiationFixedTMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvc.H"
#include "radiationModel.H"
#include "radiationConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Trad_(p.size()),
    emissivity_(0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const MarshakRadiationFixedTMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Trad_(ptf.Trad_, mapper),
    emissivity_(ptf.emissivity_)
{}


Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Trad_("Trad", dict, p.size()),
    emissivity_(readScalar(dict.lookup("emissivity")))
{
    refValue() = 4.0*radiation::sigmaSB.value()*pow4(Trad_);
    refGrad() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(refValue());
    }
}


Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const MarshakRadiationFixedTMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    Trad_(ptf.Trad_),
    emissivity_(ptf.emissivity_)
{}


Foam::MarshakRadiationFixedTMixedFvPatchScalarField::
MarshakRadiationFixedTMixedFvPatchScalarField
(
    const MarshakRadiationFixedTMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    Trad_(ptf.Trad_),
    emissivity_(ptf.emissivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MarshakRadiationFixedTMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    Trad_.autoMap(m);
}


void Foam::MarshakRadiationFixedTMixedFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const MarshakRadiationFixedTMixedFvPatchScalarField& mrptf =
        refCast<const MarshakRadiationFixedTMixedFvPatchScalarField>(ptf);

    Trad_.rmap(mrptf.Trad_, addr);
}


void Foam::MarshakRadiationFixedTMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Re-calc reference value
    refValue() = 4.0*radiation::sigmaSB.value()*pow4(Trad_);

    // Diffusion coefficient - created by radiation model's ::updateCoeffs()
    const scalarField& gamma =
        patch().lookupPatchField<volScalarField, scalar>("gammaRad");

    const scalar Ep = emissivity_/(2.0*(2.0 - emissivity_));

    // Set value fraction
    valueFraction() = 1.0/(1.0 + gamma*patch().deltaCoeffs()/Ep);

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::MarshakRadiationFixedTMixedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    Trad_.writeEntry("Trad", os);
    os.writeKeyword("emissivity") << emissivity_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        MarshakRadiationFixedTMixedFvPatchScalarField
    );
}


// ************************************************************************* //
