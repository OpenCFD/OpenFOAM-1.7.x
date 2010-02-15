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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    KinematicParcel<ParcelType>::constantProperties(parentDict),
    T0_(dimensionedScalar(this->dict().lookup("T0")).value()),
    TMin_(dimensionedScalar(this->dict().lookup("TMin")).value()),
    cp0_(dimensionedScalar(this->dict().lookup("cp0")).value()),
    epsilon0_(dimensionedScalar(this->dict().lookup("epsilon0")).value()),
    f0_(dimensionedScalar(this->dict().lookup("f0")).value()),
    Pr_(dimensionedScalar(this->dict().lookup("Pr")).value())
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::trackData::trackData
(
    ThermoCloud<ParcelType>& cloud,
    const constantProperties& constProps,
    const interpolation<scalar>& rhoInterp,
    const interpolation<vector>& UInterp,
    const interpolation<scalar>& muInterp,
    const interpolation<scalar>& TInterp,
    const interpolation<scalar>& cpInterp,
    const vector& g
)
:
    KinematicParcel<ParcelType>::trackData
    (
        cloud,
        constProps,
        rhoInterp,
        UInterp,
        muInterp,
        g
    ),
    cloud_(cloud),
    constProps_(constProps),
    TInterp_(TInterp),
    cpInterp_(cpInterp)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    ThermoCloud<ParcelType>& owner,
    const vector& position,
    const label cellI
)
:
    KinematicParcel<ParcelType>(owner, position, cellI),
    T_(0.0),
    cp_(0.0),
    Tc_(0.0),
    cpc_(0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    ThermoCloud<ParcelType>& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const constantProperties& constProps
)
:
    KinematicParcel<ParcelType>
    (
        owner,
        position,
        cellI,
        typeId,
        nParticle0,
        d0,
        U0,
        constProps
    ),
    T_(constProps.T0()),
    cp_(constProps.cp0()),
    Tc_(0.0),
    cpc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template <class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::T0() const
{
    return T0_;
}


template <class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::TMin() const
{
    return TMin_;
}


template <class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::cp0() const
{
    return cp0_;
}


template <class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::epsilon0() const
{
    return epsilon0_;
}


template <class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::f0() const
{
    return f0_;
}


template <class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::Pr() const
{
    return Pr_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ThermoCloud<ParcelType>&
Foam::ThermoParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


template <class ParcelType>
inline const typename Foam::ThermoParcel<ParcelType>::constantProperties&
Foam::ThermoParcel<ParcelType>::trackData::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackData::TInterp() const
{
    return TInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackData::cpInterp() const
{
    return cpInterp_;
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::cp() const
{
    return cp_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::T()
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::cp()
{
    return cp_;
}


// ************************************************************************* //
