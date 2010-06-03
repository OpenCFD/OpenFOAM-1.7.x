/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "DsmcParcel.H"
#include "meshTools.H"

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

    const scalar deltaT = mesh.time().deltaTValue();
    scalar tEnd = (1.0 - p.stepFraction())*deltaT;
    const scalar dtMax = tEnd;

    // For reduced-D cases, the velocity used to track needs to be
    // constrained, but the actual U_ of the parcel must not be
    // altered or used, as it is altered by patch interactions an
    // needs to retain its 3D value for collision purposes.
    vector Utracking = U_;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        // Apply correction to position for reduced-D cases
        meshTools::constrainToMeshCentre(mesh, p.position());

        Utracking = U_;

        // Apply correction to velocity to constrain tracking for
        // reduced-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);

        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        dt *= p.trackToFace(p.position() + dt*Utracking, td);

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
    label wppIndex = wpp.index();

    label wppLocalFace = wpp.whichFace(this->face());

    const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

    const scalar deltaT = td.cloud().pMesh().time().deltaTValue();

    const constantProperties& constProps(td.cloud().constProps(typeId_));

    scalar m = constProps.mass();

    vector nw = wpp.faceAreas()[wppLocalFace];
    nw /= mag(nw);

    scalar U_dot_nw = U_ & nw;

    vector Ut = U_ - U_dot_nw*nw;

    scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);

    td.cloud().rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

    td.cloud().rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

    td.cloud().linearKEBF()[wppIndex][wppLocalFace] +=
        0.5*m*(U_ & U_)*invMagUnfA;

    td.cloud().internalEBF()[wppIndex][wppLocalFace] += Ei_*invMagUnfA;

    td.cloud().iDofBF()[wppIndex][wppLocalFace] +=
        constProps.internalDegreesOfFreedom()*invMagUnfA;

    td.cloud().momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

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

    U_dot_nw = U_ & nw;

    Ut = U_ - U_dot_nw*nw;

    invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);

    td.cloud().rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

    td.cloud().rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

    td.cloud().linearKEBF()[wppIndex][wppLocalFace] +=
        0.5*m*(U_ & U_)*invMagUnfA;

    td.cloud().internalEBF()[wppIndex][wppLocalFace] += Ei_*invMagUnfA;

    td.cloud().iDofBF()[wppIndex][wppLocalFace] +=
        constProps.internalDegreesOfFreedom()*invMagUnfA;

    td.cloud().momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

    // post-interaction energy
    scalar postIE = 0.5*m*(U_ & U_) + Ei_;

    // post-interaction momentum
    vector postIMom = m*U_;

    scalar deltaQ = td.cloud().nParticle()*(preIE - postIE)/(deltaT*fA);

    vector deltaFD = td.cloud().nParticle()*(preIMom - postIMom)/(deltaT*fA);

    td.cloud().qBF()[wppIndex][wppLocalFace] += deltaQ;

    td.cloud().fDBF()[wppIndex][wppLocalFace] += deltaFD;

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
