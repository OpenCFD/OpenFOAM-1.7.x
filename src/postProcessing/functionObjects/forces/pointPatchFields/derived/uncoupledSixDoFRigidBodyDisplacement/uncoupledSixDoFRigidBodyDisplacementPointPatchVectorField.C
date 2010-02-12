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

#include "uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::
uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    motion_(),
    p0_(p.localPoints()),
    rhoInf_(1.0)
{}


uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::
uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    motion_(dict),
    rhoInf_(readScalar(dict.lookup("rhoInf")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::
uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    motion_(ptf.motion_),
    p0_(ptf.p0_, mapper),
    rhoInf_(ptf.rhoInf_)
{}


uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::
uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    motion_(ptf.motion_),
    p0_(ptf.p0_),
    rhoInf_(ptf.rhoInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();

    motion_.updatePosition(t.deltaTValue());

    vector gravity = vector::zero;

    if (db().foundObject<uniformDimensionedVectorField>("g"))
    {
        uniformDimensionedVectorField g =
        db().lookupObject<uniformDimensionedVectorField>("g");

        gravity = g.value();
    }

    // Do not modify the accelerations, except with gravity, where available
    motion_.updateForce(gravity*motion_.mass(), vector::zero, t.deltaTValue());

    Field<vector>::operator=(motion_.currentPosition(p0_) - p0_);

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    motion_.write(os);
    os.writeKeyword("rhoInf")
        << rhoInf_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
