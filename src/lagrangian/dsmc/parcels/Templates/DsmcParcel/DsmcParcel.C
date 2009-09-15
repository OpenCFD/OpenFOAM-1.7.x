/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "DsmcParcel.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::DsmcParcel<ParcelType>::move
(
    TrackData& td
)
{
    ParcelType& p = static_cast<ParcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    const scalar deltaT = td.cloud().cachedDeltaT();
    scalar tEnd = (1.0 - p.stepFraction())*deltaT;
    const scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        dt *= p.trackToFace(p.position() + dt*U_, td);

        tEnd -= dt;

        p.stepFraction() = 1.0 - tEnd/deltaT;

        if (p.onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
            {
                td.switchProcessor = true;
            }
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackData>
bool Foam::DsmcParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td,
    const label patchI
)
{
    return false;
}


template<class ParcelType>
template<class TrackData>
void Foam::DsmcParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
void Foam::DsmcParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


template<class ParcelType>
template<class TrackData>
void Foam::DsmcParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td
)
{
    const constantProperties& constProps(td.cloud().constProps(typeId_));

    scalar m = constProps.mass();

    // pre-interaction energy
    scalar preIE = 0.5*m*(U_ & U_) + Ei_;

    // pre-interaction momentum
    vector preIMom = m*U_;

    td.cloud().wallInteraction().correct
    (
        wpp,
        this->face(),
        U_,
        Ei_,
        typeId_
    );

    // post-interaction energy
    scalar postIE = 0.5*m*(U_ & U_) + Ei_;

    // post-interaction momentum
    vector postIMom = m*U_;

    label wppIndex = wpp.index();

    label wppLocalFace = wpp.whichFace(this->face());

    const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

    const scalar deltaT = td.cloud().cachedDeltaT();

    scalar deltaQ = td.cloud().nParticle()*(preIE - postIE)/(deltaT*fA);

    vector deltaFD = td.cloud().nParticle()*(preIMom - postIMom)/(deltaT*fA);

    td.cloud().q().boundaryField()[wppIndex][wppLocalFace] += deltaQ;

    td.cloud().fD().boundaryField()[wppIndex][wppLocalFace] += deltaFD;
}


template<class ParcelType>
void Foam::DsmcParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch&,
    int&
)
{}


template<class ParcelType>
template<class TrackData>
void Foam::DsmcParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParcelType>
void Foam::DsmcParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    int&
)
{}


template<class ParcelType>
void Foam::DsmcParcel<ParcelType>::transformProperties
(
    const tensor& T
)
{
    Particle<ParcelType>::transformProperties(T);
    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::DsmcParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    Particle<ParcelType>::transformProperties(separation);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "DsmcParcelIO.C"

// ************************************************************************* //

