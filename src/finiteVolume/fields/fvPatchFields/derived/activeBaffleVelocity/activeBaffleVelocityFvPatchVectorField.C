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

#include "activeBaffleVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activeBaffleVelocityFvPatchVectorField::
activeBaffleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    pName_("p"),
    cyclicPatchName_(),
    cyclicPatchLabel_(-1),
    orientation_(1),
    initWallSf_(0),
    initCyclicSf_(0),
    openFraction_(0),
    openingTime_(0),
    maxOpenFractionDelta_(0),
    curTimeIndex_(-1)
{}


Foam::activeBaffleVelocityFvPatchVectorField::
activeBaffleVelocityFvPatchVectorField
(
    const activeBaffleVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    orientation_(ptf.orientation_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    openFraction_(ptf.openFraction_),
    openingTime_(ptf.openingTime_),
    maxOpenFractionDelta_(ptf.maxOpenFractionDelta_),
    curTimeIndex_(-1)
{}


Foam::activeBaffleVelocityFvPatchVectorField::
activeBaffleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    pName_("p"),
    cyclicPatchName_(dict.lookup("cyclicPatch")),
    cyclicPatchLabel_(p.patch().boundaryMesh().findPatchID(cyclicPatchName_)),
    orientation_(readLabel(dict.lookup("orientation"))),
    initWallSf_(p.Sf()),
    initCyclicSf_(p.boundaryMesh()[cyclicPatchLabel_].Sf()),
    openFraction_(readScalar(dict.lookup("openFraction"))),
    openingTime_(readScalar(dict.lookup("openingTime"))),
    maxOpenFractionDelta_(readScalar(dict.lookup("maxOpenFractionDelta"))),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(vector::zero);

    if (dict.found("p"))
    {
        dict.lookup("p") >> pName_;
    }
}


Foam::activeBaffleVelocityFvPatchVectorField::
activeBaffleVelocityFvPatchVectorField
(
    const activeBaffleVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    orientation_(ptf.orientation_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    openFraction_(ptf.openFraction_),
    openingTime_(ptf.openingTime_),
    maxOpenFractionDelta_(ptf.maxOpenFractionDelta_),
    curTimeIndex_(-1)
{}


Foam::activeBaffleVelocityFvPatchVectorField::
activeBaffleVelocityFvPatchVectorField
(
    const activeBaffleVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    orientation_(ptf.orientation_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    openFraction_(ptf.openFraction_),
    openingTime_(ptf.openingTime_),
    maxOpenFractionDelta_(ptf.maxOpenFractionDelta_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::activeBaffleVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    //- Note: cannot map field from cyclic patch anyway so just recalculate
    //  Areas should be consistent when doing autoMap except in case of
    //  topo changes.
    //- Note: we don't want to use Sf here since triggers rebuilding of
    //  fvMesh::S() which will give problems when mapped (since already
    //  on new mesh)
    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
}

void Foam::activeBaffleVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    // See autoMap.
    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
}


void Foam::activeBaffleVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Execute the change to the openFraction only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const volScalarField& p = db().lookupObject<volScalarField>
        (
            pName_
        );

        const fvPatch& cyclicPatch = patch().boundaryMesh()[cyclicPatchLabel_];
        const labelList& cyclicFaceCells = cyclicPatch.patch().faceCells();
        label nCyclicFaces = cyclicFaceCells.size();
        label nCyclicFacesPerSide = nCyclicFaces/2;

        scalar forceDiff = 0;

        for (label facei=0; facei<nCyclicFacesPerSide; facei++)
        {
            forceDiff += p[cyclicFaceCells[facei]]*mag(initCyclicSf_[facei]);
        }

        for (label facei=nCyclicFacesPerSide; facei<nCyclicFaces; facei++)
        {
            forceDiff -= p[cyclicFaceCells[facei]]*mag(initCyclicSf_[facei]);
        }

        openFraction_ =
            max(min(
                openFraction_
              + max
                (
                    this->db().time().deltaT().value()/openingTime_,
                    maxOpenFractionDelta_
                )
               *(orientation_*sign(forceDiff)),
              1 - 1e-6), 1e-6);

        Info<< "openFraction = " << openFraction_ << endl;

        vectorField::subField Sfw = patch().patch().faceAreas();
        vectorField newSfw = (1 - openFraction_)*initWallSf_;
        forAll(Sfw, facei)
        {
            Sfw[facei] = newSfw[facei];
        }
        const_cast<scalarField&>(patch().magSf()) = mag(patch().Sf());

        const_cast<vectorField&>(cyclicPatch.Sf()) =
            openFraction_*initCyclicSf_;
        const_cast<scalarField&>(cyclicPatch.magSf()) =
            mag(cyclicPatch.Sf());

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::activeBaffleVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("cyclicPatch")
        << cyclicPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("orientation")
        << orientation_ << token::END_STATEMENT << nl;
    os.writeKeyword("openingTime")
        << openingTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("maxOpenFractionDelta")
        << maxOpenFractionDelta_ << token::END_STATEMENT << nl;
    os.writeKeyword("openFraction")
        << openFraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("p")
        << pName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        activeBaffleVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
