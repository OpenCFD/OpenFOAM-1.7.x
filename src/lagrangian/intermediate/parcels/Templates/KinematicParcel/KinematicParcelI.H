/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
inline Foam::KinematicParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    dict_(parentDict.subDict("constantProperties")),
    rhoMin_(dimensionedScalar(dict_.lookup("rhoMin")).value()),
    rho0_(dimensionedScalar(dict_.lookup("rho0")).value()),
    minParticleMass_
    (
        dimensionedScalar(dict_.lookup("minParticleMass")).value()
    )
{}


template <class ParcelType>
inline Foam::KinematicParcel<ParcelType>::trackData::trackData
(
    KinematicCloud<ParcelType>& cloud,
    const constantProperties& constProps,
    const interpolation<scalar>& rhoInterp,
    const interpolation<vector>& UInterp,
    const interpolation<scalar>& muInterp,
    const vector& g
)
:
    Particle<ParcelType>::trackData(cloud),
    cloud_(cloud),
    constProps_(constProps),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    muInterp_(muInterp),
    g_(g)
{}


template <class ParcelType>
inline Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    KinematicCloud<ParcelType>& owner,
    const vector& position,
    const label cellI
)
:
    Particle<ParcelType>(owner, position, cellI),
    typeId_(owner.parcelTypeId()),
    nParticle_(0),
    d_(0.0),
    U_(vector::zero),
    rho_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{}


template <class ParcelType>
inline Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    KinematicCloud<ParcelType>& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const constantProperties& constProps
)
:
    Particle<ParcelType>(owner, position, cellI),
    typeId_(typeId),
    nParticle_(nParticle0),
    d_(d0),
    U_(U0),
    rho_(constProps.rho0()),
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template <class ParcelType>
inline const Foam::dictionary&
Foam::KinematicParcel<ParcelType>::constantProperties::dict() const
{
    return dict_;
}


template <class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::rhoMin() const
{
    return rhoMin_;
}


template <class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_;
}


template <class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::constantProperties::minParticleMass() const
{
    return minParticleMass_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template <class ParcelType>
inline Foam::KinematicCloud<ParcelType>&
Foam::KinematicParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


template <class ParcelType>
inline const typename Foam::KinematicParcel<ParcelType>::constantProperties&
Foam::KinematicParcel<ParcelType>::trackData::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::trackData::rhoInterp() const
{
    return rhoInterp_;
}


template <class ParcelType>
inline const Foam::interpolation<Foam::vector>&
Foam::KinematicParcel<ParcelType>::trackData::UInterp() const
{
    return UInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::trackData::muInterp() const
{
    return muInterp_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::trackData::g() const
{
    return g_;
}


// * * * * * * * * * * KinematicParcel Member Functions  * * * * * * * * * * //

template <class ParcelType>
inline Foam::label Foam::KinematicParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::d() const
{
    return d_;
}


template <class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::U() const
{
    return U_;
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::rho() const
{
    return rho_;
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template <class ParcelType>
inline const Foam::vector& Foam::KinematicParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}


template <class ParcelType>
inline Foam::label Foam::KinematicParcel<ParcelType>::typeId()
{
    return typeId_;
}


template <class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::nParticle()
{
    return nParticle_;
}


template <class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::d()
{
    return d_;
}


template <class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::U()
{
    return U_;
}


template <class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::rho()
{
    return rho_;
}


template <class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::tTurb()
{
    return tTurb_;
}


template <class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::UTurb()
{
    return UTurb_;
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::wallImpactDistance
(
    const vector&
) const
{
    return 0.5*d_;
}


template <class ParcelType>
inline Foam::label Foam::KinematicParcel<ParcelType>::faceInterpolation() const
{
    // Use volume-based interpolation if dealing with external faces
    if (this->cloud().internalFace(this->face()))
    {
        return this->face();
    }
    else
    {
        return -1;
    }
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::massCell
(
    const label cellI
) const
{
    return rhoc_*this->cloud().pMesh().cellVolumes()[cellI];
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::volume() const
{
    return volume(d_);
}


template <class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::volume(const scalar d) const
{
    return mathematicalConstant::pi/6.0*pow3(d);
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaP() const
{
    return areaP(d_);
}


template <class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::areaP(const scalar d) const
{
    return 0.25*areaS(d);
}


template <class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::areaS() const
{
    return areaS(d_);
}


template <class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::areaS(const scalar d) const
{
    return mathematicalConstant::pi*d*d;
}


template <class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::Re
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar muc
) const
{
    return rhoc*mag(U - Uc_)*d/muc;
}


// ************************************************************************* //
