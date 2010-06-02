/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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

#include "timeVaryingFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingFlowRateInletVelocityFvPatchVectorField::
timeVaryingFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    flowRateInletVelocityFvPatchVectorField(p, iF),
    timeSeries_()
{}


Foam::timeVaryingFlowRateInletVelocityFvPatchVectorField::
timeVaryingFlowRateInletVelocityFvPatchVectorField
(
    const timeVaryingFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    flowRateInletVelocityFvPatchVectorField(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_)
{}


Foam::timeVaryingFlowRateInletVelocityFvPatchVectorField::
timeVaryingFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    flowRateInletVelocityFvPatchVectorField(p, iF, dict),
    timeSeries_(dict)
{}


Foam::timeVaryingFlowRateInletVelocityFvPatchVectorField::
timeVaryingFlowRateInletVelocityFvPatchVectorField
(
    const timeVaryingFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    flowRateInletVelocityFvPatchVectorField(ptf),
    timeSeries_(ptf.timeSeries_)
{}


Foam::timeVaryingFlowRateInletVelocityFvPatchVectorField::
timeVaryingFlowRateInletVelocityFvPatchVectorField
(
    const timeVaryingFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    flowRateInletVelocityFvPatchVectorField(ptf, iF),
    timeSeries_(ptf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingFlowRateInletVelocityFvPatchVectorField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    flowRate() = timeSeries_(this->db().time().timeOutputValue());
    flowRateInletVelocityFvPatchVectorField::updateCoeffs();
}


void Foam::timeVaryingFlowRateInletVelocityFvPatchVectorField::
write(Ostream& os) const
{
    flowRateInletVelocityFvPatchVectorField::write(os);
    timeSeries_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       timeVaryingFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
