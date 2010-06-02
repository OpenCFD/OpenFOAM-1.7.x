/*---------------------------------------------------------------------------* \
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

#include "SRFVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "SRFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    relative_(0),
    inletValue_(p.size(), vector::zero)
{}


SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const SRFVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    relative_(ptf.relative_),
    inletValue_(ptf.inletValue_, mapper)
{}


SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    relative_(dict.lookup("relative")),
    inletValue_("inletValue", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const SRFVelocityFvPatchVectorField& srfvpvf
)
:
    fixedValueFvPatchVectorField(srfvpvf),
    relative_(srfvpvf.relative_),
    inletValue_(srfvpvf.inletValue_)
{}


SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const SRFVelocityFvPatchVectorField& srfvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(srfvpvf, iF),
    relative_(srfvpvf.relative_),
    inletValue_(srfvpvf.inletValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SRFVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
    inletValue_.autoMap(m);
}


void SRFVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const SRFVelocityFvPatchVectorField& tiptf =
        refCast<const SRFVelocityFvPatchVectorField>(ptf);

    inletValue_.rmap(tiptf.inletValue_, addr);
}


void SRFVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If relative, include the effect of the SRF
    if (relative_)
    {
        // Get reference to the SRF model
        const SRF::SRFModel& srf =
            db().lookupObject<SRF::SRFModel>("SRFProperties");

        // Determine patch velocity due to SRF
        const vectorField SRFVelocity = srf.velocity(patch().Cf());

        operator==(-SRFVelocity + inletValue_);
    }
    // If absolute, simply supply the inlet value as a fixed value
    else
    {
        operator==(inletValue_);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void SRFVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("relative") << relative_ << token::END_STATEMENT << nl;
    inletValue_.writeEntry("inletValue", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    SRFVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
