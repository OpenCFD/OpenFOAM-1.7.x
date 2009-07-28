/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "CompositionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    mcCarrierThermo_(owner.mcCarrierThermo()),
    liquids_
    (
        liquidMixture::New
        (
            owner.mesh().objectRegistry::lookupObject<dictionary>
            (
                owner.carrierThermo().name()
            )
        )
    ),
    solids_
    (
        solidMixture::New
        (
            owner.mesh().objectRegistry::lookupObject<dictionary>
            (
                owner.carrierThermo().name()
            )
        )
    ),
    phaseProps_
    (
        coeffDict_.lookup("phases"),
        mcCarrierThermo_.species(),
        liquids_().components(),
        solids_().components()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::~CompositionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::CompositionModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::CompositionModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary& Foam::CompositionModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
const Foam::multiComponentMixture<typename CloudType::thermoType>&
Foam::CompositionModel<CloudType>::mcCarrierThermo() const
{
    return mcCarrierThermo_;
}


template<class CloudType>
const Foam::liquidMixture& Foam::CompositionModel<CloudType>::liquids() const
{
    return liquids_();
}


template<class CloudType>
const Foam::solidMixture& Foam::CompositionModel<CloudType>::solids() const
{
    return solids_();
}


template<class CloudType>
const Foam::phasePropertiesList&
Foam::CompositionModel<CloudType>::phaseProps() const
{
    return phaseProps_;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::nPhase() const
{
    return phaseProps_.size();
}


template<class CloudType>
const Foam::wordList&
Foam::CompositionModel<CloudType>::phaseTypes() const
{
    // if only 1 phase, return the constituent component names
    if (phaseProps_.size() == 1)
    {
        return phaseProps_[0].names();
    }
    else
    {
        return phaseProps_.phaseTypes();
    }
}


template<class CloudType>
const Foam::wordList&
Foam::CompositionModel<CloudType>::stateLabels() const
{
    return phaseProps_.stateLabels();
}


template<class CloudType>
const Foam::wordList&
Foam::CompositionModel<CloudType>::componentNames(const label phaseI) const
{
    return phaseProps_[phaseI].names();
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::globalCarrierId
(
    const word& cmptName
) const
{
    forAll(mcCarrierThermo_.species(), i)
    {
        if (cmptName == mcCarrierThermo_.species()[i])
        {
            return i;
        }
    }

    FatalErrorIn
    (
        "Foam::label Foam::CompositionModel<CloudType>::globalCarrierId"
        "("
            "const word&"
        ") const"
    )   << "Unable to determine global id for requested component "
        << cmptName << nl << abort(FatalError);

    return -1;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::globalId
(
    const label phaseI,
    const word& cmptName
) const
{
    label id = phaseProps_[phaseI].globalId(cmptName);

    if (id < 0)
    {
        FatalErrorIn
        (
            "Foam::label Foam::CompositionModel<CloudType>::globalId"
            "("
                "const label, "
                "const word&"
            ") const"
        )   << "Unable to determine global id for requested component "
            << cmptName << nl << abort(FatalError);
    }

    return id;
}


template<class CloudType>
const Foam::labelList& Foam::CompositionModel<CloudType>::globalIds
(
    const label phaseI
) const
{
    return phaseProps_[phaseI].globalIds();
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::localId
(
    const label phaseI,
    const word& cmptName
) const
{
    label id = phaseProps_[phaseI].id(cmptName);

    if (id < 0)
    {
        FatalErrorIn
        (
            "Foam::label Foam::CompositionModel<CloudType>::localId"
            "("
                "const label, "
                "const word&"
            ") const"
        )   << "Unable to determine local id for component " << cmptName
            << nl << abort(FatalError);
    }

    return id;
}


template<class CloudType>
Foam::label Foam::CompositionModel<CloudType>::localToGlobalCarrierId
(
    const label phaseI,
    const label id
) const
{
    label gid = phaseProps_[phaseI].globalCarrierIds()[id];

    if (gid < 0)
    {
        FatalErrorIn
        (
            "Foam::label "
            "Foam::CompositionModel<CloudType>::localToGlobalCarrierId"
            "("
                "const label, "
                "const label"
            ") const"
        )   << "Unable to determine global carrier id for phase "
            << phaseI << " with local id " << id
            << nl << abort(FatalError);
    }

    return gid;
}


template<class CloudType>
const Foam::scalarField& Foam::CompositionModel<CloudType>::Y0
(
    const label phaseI
) const
{
    return phaseProps_[phaseI].Y();
}


template<class CloudType>
Foam::scalarField Foam::CompositionModel<CloudType>::X
(
    const label phaseI,
    const scalarField& Y
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalarField X(Y.size(), 0.0);
    scalar WInv = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                WInv += Y[i]/mcCarrierThermo_.speciesData()[gid].W();
                X[i] = Y[i]/mcCarrierThermo_.speciesData()[gid].W();
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                WInv += Y[i]/this->liquids().properties()[gid].W();
                X[i] += Y[i]/this->liquids().properties()[gid].W();
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalarField Foam::CompositionModel<CloudType>::X"
                "("
                    "const label, "
                    "const scalarField&"
                ") const"
            )   << "Only possible to convert gas and liquid mass fractions"
                << nl << abort(FatalError);
        }
    }

    return X/WInv;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::H
(
    const label phaseI,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalar HMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HMixture += Y[i]*mcCarrierThermo_.speciesData()[gid].H(T);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HMixture += Y[i]*this->liquids().properties()[gid].h(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                HMixture +=
                     Y[i]
                    *(
                        this->solids().properties()[gid].Hf()
                      + this->solids().properties()[gid].cp()*T
                     );
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::CompositionModel<CloudType>::H"
                "("
                "    const label, "
                "    const scalarField&, "
                "    const scalar, "
                "    const scalar"
                ") const"
            )   << "Unknown phase enumeration" << nl << abort(FatalError);
        }
    }

    return HMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::cp
(
    const label phaseI,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalar cpMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                cpMixture += Y[i]*mcCarrierThermo_.speciesData()[gid].Cp(T);
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                cpMixture += Y[i]*this->liquids().properties()[gid].cp(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                cpMixture += Y[i]*this->solids().properties()[gid].cp();
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::CompositionModel<CloudType>::cp"
                "("
                    "const label, "
                    "const scalarField&, "
                    "const scalar, "
                    "const scalar"
                ") const"
            )   << "Unknown phase enumeration" << nl << abort(FatalError);
        }
    }

    return cpMixture;
}


template<class CloudType>
Foam::scalar Foam::CompositionModel<CloudType>::L
(
    const label phaseI,
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    const phaseProperties& props = phaseProps_[phaseI];
    scalar LMixture = 0.0;
    switch (props.phase())
    {
        case phaseProperties::GAS:
        {
            if (debug)
            {
                WarningIn
                (
                    "Foam::scalar Foam::CompositionModel<CloudType>::L"
                    "("
                        "const label, "
                        "const scalarField&, "
                        "const scalar, "
                        "const scalar"
                    ") const\n"
                )   << "No support for gaseous components" << endl;
            }
            break;
        }
        case phaseProperties::LIQUID:
        {
            forAll(Y, i)
            {
                label gid = props.globalIds()[i];
                LMixture += Y[i]*this->liquids().properties()[gid].hl(p, T);
            }
            break;
        }
        case phaseProperties::SOLID:
        {
            if (debug)
            {
                WarningIn
                (
                    "Foam::scalar Foam::CompositionModel<CloudType>::L"
                    "("
                        "const label, "
                        "const scalarField&, "
                        "const scalar, "
                        "const scalar"
                    ") const\n"
                )   << "No support for solid components" << endl;
            }
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::CompositionModel<CloudType>::L"
                "("
                    "const label, "
                    "const scalarField&, "
                    "const scalar, "
                    "const scalar"
                ") const"
            )   << "Unknown phase enumeration" << nl << abort(FatalError);
        }
    }

    return LMixture;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewCompositionModel.C"

// ************************************************************************* //

