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

#include "cylindricalInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    axialVelocity_(0),
    centre_(pTraits<vector>::zero),
    axis_(pTraits<vector>::zero),
    rpm_(0),
    radialVelocity_(0)
{}


Foam::
cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    axialVelocity_(ptf.axialVelocity_),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    rpm_(ptf.rpm_),
    radialVelocity_(ptf.radialVelocity_)
{}


Foam::
cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    axialVelocity_(readScalar(dict.lookup("axialVelocity"))),
    centre_(dict.lookup("centre")),
    axis_(dict.lookup("axis")),
    rpm_(readScalar(dict.lookup("rpm"))),
    radialVelocity_(readScalar(dict.lookup("radialVelocity")))
{}


Foam::
cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    axialVelocity_(ptf.axialVelocity_),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    rpm_(ptf.rpm_),
    radialVelocity_(ptf.radialVelocity_)
{}


Foam::
cylindricalInletVelocityFvPatchVectorField::
cylindricalInletVelocityFvPatchVectorField
(
    const cylindricalInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    axialVelocity_(ptf.axialVelocity_),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    rpm_(ptf.rpm_),
    radialVelocity_(ptf.radialVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylindricalInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vector hatAxis = axis_/mag(axis_);

    vectorField r = (patch().Cf() - centre_);

    vectorField d =  r - (hatAxis & r)*hatAxis;

    vectorField tangVelo =
        (rpm_*mathematicalConstant::pi/30.0)*(hatAxis)^d;

    operator==(tangVelo + axis_*axialVelocity_ + radialVelocity_*d);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::cylindricalInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("axialVelocity") << axialVelocity_ <<
        token::END_STATEMENT << nl;
    os.writeKeyword("centre") << centre_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("rpm") << rpm_ << token::END_STATEMENT << nl;
    os.writeKeyword("radialVelocity") << radialVelocity_ <<
        token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       cylindricalInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
