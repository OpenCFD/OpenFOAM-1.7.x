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

#include "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "directMappedPatchBase.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::interfaceOwner
(
    const polyMesh& nbrRegion,
    const polyPatch& nbrPatch
) const
{
    const fvMesh& myRegion = patch().boundaryMesh().mesh();

    if (nbrRegion.name() == myRegion.name())
    {
        return patch().index() < nbrPatch.index();
    }
    else
    {
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
                    "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField"
                    "::interfaceOwner(const polyMesh&"
                    ", const polyPatch&)const"
                )   << "Cannot find region " << myRegion.name()
                    << " neither in fluids " << props.fluidRegionNames()
                    << " nor in solids " << props.solidRegionNames()
                    << exit(FatalError);
            }
            myIndex = props.fluidRegionNames().size() + i;
        }
        label nbrIndex = findIndex
        (
            props.fluidRegionNames(),
            nbrRegion.name()
        );
        if (nbrIndex == -1)
        {
            label i = findIndex(props.solidRegionNames(), nbrRegion.name());

            if (i == -1)
            {
                FatalErrorIn
                (
                    "coupleManager::interfaceOwner"
                    "(const polyMesh&, const polyPatch&) const"
                )   << "Cannot find region " << nbrRegion.name()
                    << " neither in fluids " << props.fluidRegionNames()
                    << " nor in solids " << props.solidRegionNames()
                    << exit(FatalError);
            }
            nbrIndex = props.fluidRegionNames().size() + i;
        }
        return myIndex < nbrIndex;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
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


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& ptf,
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


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
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
            "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::"
            "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField\n"
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


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    KName_(wtcsf.KName_),
    fixesValue_(wtcsf.fixesValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::K() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (KName_ == "none")
    {
        const compressible::RASModel& model =
            db().lookupObject<compressible::RASModel>("RASProperties");

        tmp<volScalarField> talpha = model.alphaEff();

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        return
            talpha().boundaryField()[patch().index()]
           *thermo.Cp()().boundaryField()[patch().index()];
    }
    else if (mesh.objectRegistry::foundObject<volScalarField>(KName_))
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
            "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::K()"
            " const"
        )   << "Did not find field " << KName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'K' to 'none', a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::updateCoeffs()
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
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    // Force recalculation of mapping and schedule
    const mapDistribute& distMap = mpp.map();
    (void)distMap.schedule();

    tmp<scalarField> intFld = patchInternalField();

    if (interfaceOwner(nbrMesh, nbrPatch.patch()))
    {
        // Note: other side information could be cached - it only needs
        // to be updated the first time round the iteration (i.e. when
        // switching regions) but unfortunately we don't have this information.


        // Calculate the temperature by harmonic averaging
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField&
        nbrField =
        refCast
        <
            const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
        >
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
        const_cast<turbulentTemperatureCoupledBaffleMixedFvPatchScalarField&>
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

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " -> "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
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
            this->refGrad()[i] = 0.0;   // not used by me
            this->valueFraction()[i] = 1.0;
            nFixed++;
        }
        else
        {
            // Fixed gradient. Make sure to have valid refValue (even though
            // I am not using it - other boundary conditions might)
            this->refValue()[i] = operator[](i);
            this->refGrad()[i] = normalGradient()[i];
            this->valueFraction()[i] = 0.0;
        }
    }

    reduce(nFixed, sumOp<label>());

    fixesValue_ = (nFixed > 0);

    if (debug)
    {
        label nTotSize = returnReduce(this->size(), sumOp<label>());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " -> "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " patch:" << patch().name()
            << " out of:" << nTotSize
            << " fixedBC:" << nFixed
            << " gradient:" << nTotSize-nFixed << endl;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::write
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

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
