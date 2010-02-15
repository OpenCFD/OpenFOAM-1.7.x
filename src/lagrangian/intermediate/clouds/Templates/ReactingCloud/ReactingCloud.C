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

#include "ReactingCloud.H"

#include "CompositionModel.H"
#include "PhaseChangeModel.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::checkSuppliedComposition
(
    const scalarField& YSupplied,
    const scalarField& Y,
    const word& YName
)
{
    if (YSupplied.size() != Y.size())
    {
        FatalErrorIn
        (
            "ReactingCloud<ParcelType>::checkSuppliedComposition"
            "("
                "const scalarField&, "
                "const scalarField&, "
                "const word&"
            ")"
        )   << YName << " supplied, but size is not compatible with "
            << "parcel composition: " << nl << "    "
            << YName << "(" << YSupplied.size() << ") vs required composition "
            << YName << "(" << Y.size() << ")" << nl
            << abort(FatalError);
    }
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::preEvolve()
{
    ThermoCloud<ParcelType>::preEvolve();
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::evolveCloud()
{
    const volScalarField& T = this->carrierThermo().T();
    const volScalarField cp = this->carrierThermo().Cp();
    const volScalarField& p = this->carrierThermo().p();

    autoPtr<interpolation<scalar> > rhoInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->rho()
    );

    autoPtr<interpolation<vector> > UInterp = interpolation<vector>::New
    (
        this->interpolationSchemes(),
        this->U()
    );

    autoPtr<interpolation<scalar> > muInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->mu()
    );

    autoPtr<interpolation<scalar> > TInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        T
    );

    autoPtr<interpolation<scalar> > cpInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        cp
    );

    autoPtr<interpolation<scalar> > pInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        p
    );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterp(),
        UInterp(),
        muInterp(),
        TInterp(),
        cpInterp(),
        pInterp(),
        this->g().value()
    );

    this->injection().inject(td);

    if (this->coupled())
    {
        resetSourceTerms();
    }

    Cloud<ParcelType>::move(td);
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::postEvolve()
{
    ThermoCloud<ParcelType>::postEvolve();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingCloud<ParcelType>::ReactingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    basicThermo& thermo,
    bool readFields
)
:
    ThermoCloud<ParcelType>(cloudName, rho, U, g, thermo, false),
    reactingCloud(),
    constProps_(this->particleProperties()),
    mcCarrierThermo_
    (
        dynamic_cast<multiComponentMixture<thermoType>&>(thermo)
    ),
    compositionModel_
    (
        CompositionModel<ReactingCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    phaseChangeModel_
    (
        PhaseChangeModel<ReactingCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    rhoTrans_(mcCarrierThermo_.species().size()),
    dMassPhaseChange_(0.0)
{
    // Set storage for mass source fields and initialise to zero
    forAll(rhoTrans_, i)
    {
        rhoTrans_.set
        (
            i,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    this->name() + "rhoTrans_" + mcCarrierThermo_.species()[i],
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
    }

    if (readFields)
    {
        ParcelType::readFields(*this);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingCloud<ParcelType>::~ReactingCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    ThermoCloud<ParcelType>::checkParcelProperties
    (
        parcel,
        lagrangianDt,
        fullyDescribed
    );

    if (!fullyDescribed)
    {
        parcel.Y() = composition().YMixture0();
    }
    else
    {
        checkSuppliedComposition
        (
            parcel.Y(),
            composition().YMixture0(),
            "YMixture"
        );
    }

    // derived information - store initial mass
    parcel.mass0() = parcel.mass();
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::resetSourceTerms()
{
    ThermoCloud<ParcelType>::resetSourceTerms();
    forAll(rhoTrans_, i)
    {
        rhoTrans_[i].field() = 0.0;
    }
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::evolve()
{
    if (this->active())
    {
        preEvolve();

        evolveCloud();

        postEvolve();

        info();
        Info<< endl;
    }
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::info() const
{
    ThermoCloud<ParcelType>::info();

    Info<< "    Mass transfer phase change      = "
        << returnReduce(dMassPhaseChange_, sumOp<scalar>()) << nl;
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::addToMassPhaseChange(const scalar dMass)
{
    dMassPhaseChange_ += dMass;
}


// ************************************************************************* //
