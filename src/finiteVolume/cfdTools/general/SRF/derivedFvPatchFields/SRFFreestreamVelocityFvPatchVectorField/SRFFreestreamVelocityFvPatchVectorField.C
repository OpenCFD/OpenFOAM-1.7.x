/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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

#include "SRFFreestreamVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "SRFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SRFFreestreamVelocityFvPatchVectorField::SRFFreestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(p, iF),
    UInf_(vector::zero)
{}


SRFFreestreamVelocityFvPatchVectorField::SRFFreestreamVelocityFvPatchVectorField
(
    const SRFFreestreamVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchVectorField(ptf, p, iF, mapper),
    UInf_(ptf.UInf_)
{}


SRFFreestreamVelocityFvPatchVectorField::SRFFreestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchVectorField(p, iF),
    UInf_(dict.lookup("UInf"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


SRFFreestreamVelocityFvPatchVectorField::SRFFreestreamVelocityFvPatchVectorField
(
    const SRFFreestreamVelocityFvPatchVectorField& srfvpvf
)
:
    inletOutletFvPatchVectorField(srfvpvf),
    UInf_(srfvpvf.UInf_)
{}


SRFFreestreamVelocityFvPatchVectorField::SRFFreestreamVelocityFvPatchVectorField
(
    const SRFFreestreamVelocityFvPatchVectorField& srfvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(srfvpvf, iF),
    UInf_(srfvpvf.UInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SRFFreestreamVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get reference to the SRF model
    const SRF::SRFModel& srf =
        db().lookupObject<SRF::SRFModel>("SRFProperties");

    scalar time = this->db().time().value();
    scalar theta = time*mag(srf.omega().value());

    refValue() =
        cos(theta)*UInf_ + sin(theta)*(srf.axis() ^ UInf_)
      - srf.velocity(patch().Cf());

    // Set the inlet-outlet choice based on the direction of the freestream
    valueFraction() = 1.0 - pos(refValue() & patch().Sf());

    mixedFvPatchField<vector>::updateCoeffs();
}


void SRFFreestreamVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("UInf") << UInf_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    SRFFreestreamVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
