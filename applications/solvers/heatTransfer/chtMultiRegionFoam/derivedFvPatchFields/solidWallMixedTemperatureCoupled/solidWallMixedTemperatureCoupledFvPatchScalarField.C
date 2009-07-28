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

#include "solidWallMixedTemperatureCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "directMappedPatchBase.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::interfaceOwner
(
    const polyMesh& nbrRegion
) const
{
    const fvMesh& myRegion = patch().boundaryMesh().mesh();

    const regionProperties& props =
        myRegion.objectRegistry::parent().lookupObject<regionProperties>
        (
            "regionProperties"
        );

    label myIndex = findIndex(props.fluidRegionNames(), myRegion.name());
    if (myIndex == -1)
    {
        label i = findIndex(props.solidRegionNames(), myRegion.name());

        if (i == -1)
        {
            FatalErrorIn
            (
                "solidWallMixedTemperatureCoupledFvPatchScalarField"
                "::interfaceOwner(const polyMesh&) const"
            )   << "Cannot find region " << myRegion.name()
                << " neither in fluids " << props.fluidRegionNames()
                << " nor in solids " << props.solidRegionNames()
                << exit(FatalError);
        }
        myIndex = props.fluidRegionNames().size() + i;
    }
    label nbrIndex = findIndex(props.fluidRegionNames(), nbrRegion.name());
    if (nbrIndex == -1)
    {
        label i = findIndex(props.solidRegionNames(), nbrRegion.name());

        if (i == -1)
        {
            FatalErrorIn("coupleManager::interfaceOwner(const polyMesh&) const")
                << "Cannot find region " << nbrRegion.name()
                << " neither in fluids " << props.fluidRegionNames()
                << " nor in solids " << props.solidRegionNames()
                << exit(FatalError);
        }
        nbrIndex = props.fluidRegionNames().size() + i;
    }

    return myIndex < nbrIndex;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName"),
    KName_("undefined-K")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
    this->fixesValue_ = true;
}


Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const solidWallMixedTemperatureCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_),
    KName_(ptf.KName_),
    fixesValue_(ptf.fixesValue_)
{}


Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    KName_(dict.lookup("K"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "solidWallMixedTemperatureCoupledFvPatchScalarField::"
            "solidWallMixedTemperatureCoupledFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
        fixesValue_ = readBool(dict.lookup("fixesValue"));
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
        fixesValue_ = true;
    }
}


Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const solidWallMixedTemperatureCoupledFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    KName_(wtcsf.KName_),
    fixesValue_(wtcsf.fixesValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvPatchScalarField&
Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::K() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(KName_);
}


void Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the directMappedPatchBase
    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();

    tmp<scalarField> intFld = patchInternalField();

    if (interfaceOwner(nbrMesh))
    {
        // Note: other side information could be cached - it only needs
        // to be updated the first time round the iteration (i.e. when
        // switching regions) but unfortunately we don't have this information.

        const mapDistribute& distMap = mpp.map();
        const fvPatch& nbrPatch = refCast<const fvMesh>
        (
            nbrMesh
        ).boundary()[mpp.samplePolyPatch().index()];


        // Calculate the temperature by harmonic averaging
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const solidWallMixedTemperatureCoupledFvPatchScalarField& nbrField =
        refCast<const solidWallMixedTemperatureCoupledFvPatchScalarField>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldName_
            )
        );

        // Swap to obtain full local values of neighbour internal field
        scalarField nbrIntFld = nbrField.patchInternalField();
        mapDistribute::distribute
        (
            Pstream::defaultCommsType,
            distMap.schedule(),
            distMap.constructSize(),
            distMap.subMap(),           // what to send
            distMap.constructMap(),     // what to receive
            nbrIntFld
        );

        // Swap to obtain full local values of neighbour K*delta
        scalarField nbrKDelta = nbrField.K()*nbrPatch.deltaCoeffs();
        mapDistribute::distribute
        (
            Pstream::defaultCommsType,
            distMap.schedule(),
            distMap.constructSize(),
            distMap.subMap(),           // what to send
            distMap.constructMap(),     // what to receive
            nbrKDelta
        );


        tmp<scalarField> myKDelta = K()*patch().deltaCoeffs();

        // Calculate common wall temperature. Reuse *this to store common value.
        scalarField Twall
        (
            (myKDelta()*intFld() + nbrKDelta*nbrIntFld)
          / (myKDelta() + nbrKDelta)
        );
        // Assign to me
        fvPatchScalarField::operator=(Twall);
        // Distribute back and assign to neighbour
        mapDistribute::distribute
        (
            Pstream::defaultCommsType,
            distMap.schedule(),
            nbrField.size(),
            distMap.constructMap(),     // reverse : what to send
            distMap.subMap(),
            Twall
        );
        const_cast<solidWallMixedTemperatureCoupledFvPatchScalarField&>
        (
            nbrField
        ).fvPatchScalarField::operator=(Twall);
    }


    // Switch between fixed value (of harmonic avg) or gradient
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nFixed = 0;

    // Like snGrad but bypass switching on refValue/refGrad.
    tmp<scalarField> normalGradient = (*this-intFld())*patch().deltaCoeffs();

    if (debug)
    {
        scalar Q = gSum(K()*patch().magSf()*normalGradient());

        Info<< "solidWallMixedTemperatureCoupledFvPatchScalarField::"
            << "updateCoeffs() :"
            << " patch:" << patch().name()
            << " heatFlux:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    forAll(*this, i)
    {
        // if outgoing flux use fixed value.
        if (normalGradient()[i] < 0.0)
        {
            this->refValue()[i] = operator[](i);
            this->refGrad()[i] = 0.0;   // not used
            this->valueFraction()[i] = 1.0;
            nFixed++;
        }
        else
        {
            this->refValue()[i] = 0.0;  // not used
            this->refGrad()[i] = normalGradient()[i];
            this->valueFraction()[i] = 0.0;
        }
    }

    reduce(nFixed, sumOp<label>());

    fixesValue_ = (nFixed > 0);

    if (debug)
    {
        label nTotSize = returnReduce(this->size(), sumOp<label>());

        Info<< "solidWallMixedTemperatureCoupledFvPatchScalarField::"
            << "updateCoeffs() :"
            << " patch:" << patch().name()
            << " out of:" << nTotSize
            << " fixedBC:" << nFixed
            << " gradient:" << nTotSize-nFixed << endl;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fixesValue") << fixesValue_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    solidWallMixedTemperatureCoupledFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
