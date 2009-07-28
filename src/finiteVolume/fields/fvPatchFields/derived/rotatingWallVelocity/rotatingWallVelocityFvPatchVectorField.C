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

#include "rotatingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rotatingWallVelocityFvPatchVectorField::rotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(vector::zero),
    axis_(vector::zero),
    omega_(0)
{}


rotatingWallVelocityFvPatchVectorField::rotatingWallVelocityFvPatchVectorField
(
    const rotatingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_)
{}


rotatingWallVelocityFvPatchVectorField::rotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    omega_(readScalar(dict.lookup("omega")))
{
    // Evaluate the wall velocity
    updateCoeffs();
}


rotatingWallVelocityFvPatchVectorField::rotatingWallVelocityFvPatchVectorField
(
    const rotatingWallVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchField<vector>(pivpvf),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    omega_(pivpvf.omega_)
{}


rotatingWallVelocityFvPatchVectorField::rotatingWallVelocityFvPatchVectorField
(
    const rotatingWallVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(pivpvf, iF),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    omega_(pivpvf.omega_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rotatingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Calculate the rotating wall velocity from the specification of the motion
    vectorField Up = (-omega_)*((patch().Cf() - origin_) ^ (axis_/mag(axis_)));

    // Remove the component of Up normal to the wall
    // just in case it is not exactly circular
    vectorField n = patch().nf();
    vectorField::operator=(Up - n*(n & Up));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void rotatingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega") << omega_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    rotatingWallVelocityFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
