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

#include "Particle.H"
#include "Cloud.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "wallPolyPatch.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::Particle<ParticleType>::findFaces
(
    const vector& position,
    DynamicList<label>& faceList
) const
{
    const polyMesh& mesh = cloud_.polyMesh_;
    const labelList& faces = mesh.cells()[celli_];
    const vector& C = mesh.cellCentres()[celli_];

    faceList.clear();
    forAll(faces, i)
    {
        label facei = faces[i];
        scalar lam = lambda(C, position, facei);

        if ((lam > 0) && (lam < 1.0))
        {
            faceList.append(facei);
        }
    }
}


template<class ParticleType>
void Foam::Particle<ParticleType>::findFaces
(
    const vector& position,
    const label celli,
    const scalar stepFraction,
    DynamicList<label>& faceList
) const
{
    const polyMesh& mesh = cloud_.pMesh();
    const labelList& faces = mesh.cells()[celli];
    const vector& C = mesh.cellCentres()[celli];

    faceList.clear();
    forAll(faces, i)
    {
        label facei = faces[i];
        scalar lam = lambda(C, position, facei, stepFraction);

        if ((lam > 0) && (lam < 1.0))
        {
            faceList.append(facei);
        }
    }
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::prepareForParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    // Convert the face index to be local to the processor patch
    facei_ = patchFace(patchi, facei_);
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::correctAfterParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    const processorPolyPatch& ppp =
        refCast<const processorPolyPatch>
        (cloud_.pMesh().boundaryMesh()[patchi]);

    celli_ = ppp.faceCells()[facei_];

    if (!ppp.parallel())
    {
        if (ppp.forwardT().size() == 1)
        {
            const tensor& T = ppp.forwardT()[0];
            transformPosition(T);
            static_cast<ParticleType&>(*this).transformProperties(T);
        }
        else
        {
            const tensor& T = ppp.forwardT()[facei_];
            transformPosition(T);
            static_cast<ParticleType&>(*this).transformProperties(T);
        }
    }
    else if (ppp.separated())
    {
        if (ppp.separation().size() == 1)
        {
            position_ -= ppp.separation()[0];
            static_cast<ParticleType&>(*this).transformProperties
            (
                -ppp.separation()[0]
            );
        }
        else
        {
            position_ -= ppp.separation()[facei_];
            static_cast<ParticleType&>(*this).transformProperties
            (
                -ppp.separation()[facei_]
            );
        }
    }

    // Reset the face index for the next tracking operation
    if (stepFraction_ > (1.0 - SMALL))
    {
        stepFraction_ = 1.0;
        facei_ = -1;
    }
    else
    {
        facei_ += ppp.start();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Particle<ParticleType>::Particle
(
    const Cloud<ParticleType>& cloud,
    const vector& position,
    const label celli
)
:
    cloud_(cloud),
    position_(position),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    origProc_(Pstream::myProcNo()),
    origId_(cloud_.getNewParticleID())
{}


template<class ParticleType>
Foam::Particle<ParticleType>::Particle(const Particle<ParticleType>& p)
:
    cloud_(p.cloud_),
    position_(p.position_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackData>
Foam::label Foam::Particle<ParticleType>::track
(
    const vector& endPosition,
    TrackData& td
)
{
    facei_ = -1;

    // Tracks to endPosition or stop on boundary
    while(!onBoundary() && stepFraction_ < 1.0 - SMALL)
    {
        stepFraction_ += trackToFace(endPosition, td)*(1.0 - stepFraction_);
    }

    return facei_;
}



template<class ParticleType>
Foam::label Foam::Particle<ParticleType>::track(const vector& endPosition)
{
    int dummyTd;
    return track(endPosition, dummyTd);
}

template<class ParticleType>
template<class TrackData>
Foam::scalar Foam::Particle<ParticleType>::trackToFace
(
    const vector& endPosition,
    TrackData& td
)
{
    const polyMesh& mesh = cloud_.polyMesh_;

    DynamicList<label>& faces = cloud_.labels_;
    findFaces(endPosition, faces);

    facei_ = -1;
    scalar trackFraction = 0.0;

    if (faces.empty())  // inside cell
    {
        trackFraction = 1.0;
        position_ = endPosition;
    }
    else // hit face
    {
        scalar lambdaMin = GREAT;

        if (faces.size() == 1)
        {
            lambdaMin = lambda(position_, endPosition, faces[0], stepFraction_);
            facei_ = faces[0];
        }
        else
        {
            // If the particle has to cross more than one cell to reach the
            // endPosition, we check which way to go.
            // If one of the faces is a boundary face and the particle is
            // outside, we choose the boundary face.
            // The particle is outside if one of the lambda's is > 1 or < 0
            forAll(faces, i)
            {
                scalar lam =
                    lambda(position_, endPosition, faces[i], stepFraction_);

                if (lam < lambdaMin)
                {
                    lambdaMin = lam;
                    facei_ = faces[i];
                }
            }
        }

        bool internalFace = cloud_.internalFace(facei_);

        // For warped faces the particle can be 'outside' the cell.
        // This will yield a lambda larger than 1, or smaller than 0
        // For values < 0, the particle travels away from the cell
        // and we don't move the particle, only change cell.
        // For values larger than 1, we move the particle to endPosition only.
        if (lambdaMin > 0.0)
        {
            if (lambdaMin <= 1.0)
            {
                trackFraction = lambdaMin;
                position_ += trackFraction*(endPosition - position_);
            }
            else
            {
                trackFraction = 1.0;
                position_ = endPosition;
            }
        }
        else if (static_cast<ParticleType&>(*this).softImpact())
        {
            // Soft-sphere particles can travel outside the domain
            // but we don't use lambda since this the particle
            // is going away from face
            trackFraction = 1.0;
            position_ = endPosition;
        }

        // change cell
        if (internalFace) // Internal face
        {
            if (celli_ == mesh.faceOwner()[facei_])
            {
                celli_ = mesh.faceNeighbour()[facei_];
            }
            else if (celli_ == mesh.faceNeighbour()[facei_])
            {
                celli_ = mesh.faceOwner()[facei_];
            }
            else
            {
                FatalErrorIn
                (
                    "Particle::trackToFace(const vector&, TrackData&)"
                )<< "addressing failure" << nl
                 << abort(FatalError);
            }
        }
        else
        {
            ParticleType& p = static_cast<ParticleType&>(*this);

            // Soft-sphere algorithm ignores the boundary
            if (p.softImpact())
            {
                trackFraction = 1.0;
                position_ = endPosition;
            }

            label patchi = patch(facei_);
            const polyPatch& patch = mesh.boundaryMesh()[patchi];

            if (!p.hitPatch(patch, td, patchi))
            {
                if (isA<wedgePolyPatch>(patch))
                {
                    p.hitWedgePatch
                    (
                        static_cast<const wedgePolyPatch&>(patch), td
                    );
                }
                else if (isA<symmetryPolyPatch>(patch))
                {
                    p.hitSymmetryPatch
                    (
                        static_cast<const symmetryPolyPatch&>(patch), td
                    );
                }
                else if (isA<cyclicPolyPatch>(patch))
                {
                    p.hitCyclicPatch
                    (
                        static_cast<const cyclicPolyPatch&>(patch), td
                    );
                }
                else if (isA<processorPolyPatch>(patch))
                {
                    p.hitProcessorPatch
                    (
                        static_cast<const processorPolyPatch&>(patch), td
                    );
                }
                else if (isA<wallPolyPatch>(patch))
                {
                    p.hitWallPatch
                    (
                        static_cast<const wallPolyPatch&>(patch), td
                    );
                }
                else
                {
                    p.hitPatch(patch, td);
                }
            }
        }
    }

    // If the trackFraction = 0 something went wrong.
    // Either the particle is flipping back and forth across a face perhaps
    // due to velocity interpolation errors or it is in a "hole" in the mesh
    // caused by face warpage.
    // In both cases resolve the positional ambiguity by moving the particle
    // slightly towards the cell-centre.
    if (trackFraction < SMALL)
    {
        position_ += 1.0e-3*(mesh.cellCentres()[celli_] - position_);
    }

    return trackFraction;
}

template<class ParticleType>
Foam::scalar Foam::Particle<ParticleType>::trackToFace
(
    const vector& endPosition
)
{
    int dummyTd;
    return trackToFace(endPosition, dummyTd);
}

template<class ParticleType>
void Foam::Particle<ParticleType>::transformPosition(const tensor& T)
{
    position_ = transform(T, position_);
}


template<class ParticleType>
void Foam::Particle<ParticleType>::transformProperties(const tensor&)
{}


template<class ParticleType>
void Foam::Particle<ParticleType>::transformProperties(const vector&)
{}


template<class ParticleType>
template<class TrackData>
bool Foam::Particle<ParticleType>::hitPatch
(
    const polyPatch&,
    TrackData&,
    const label
)
{
    return false;
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitWedgePatch
(
    const wedgePolyPatch& wpp,
    TrackData&
)
{
    vector nf = wpp.faceAreas()[wpp.whichFace(facei_)];
    nf /= mag(nf);

    static_cast<ParticleType&>(*this).transformProperties(I - 2.0*nf*nf);
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitSymmetryPatch
(
    const symmetryPolyPatch& spp,
    TrackData&
)
{
    vector nf = spp.faceAreas()[spp.whichFace(facei_)];
    nf /= mag(nf);

    static_cast<ParticleType&>(*this).transformProperties(I - 2.0*nf*nf);
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitCyclicPatch
(
    const cyclicPolyPatch& cpp,
    TrackData&
)
{
    label patchFacei_ = cpp.whichFace(facei_);

    facei_ = cpp.transformGlobalFace(facei_);

    celli_ = cloud_.polyMesh_.faceOwner()[facei_];

    if (!cpp.parallel())
    {
        const tensor& T = cpp.transformT(patchFacei_);

        transformPosition(T);
        static_cast<ParticleType&>(*this).transformProperties(T);
    }
    else if (cpp.separated())
    {
        position_ += cpp.separation(patchFacei_);
        static_cast<ParticleType&>(*this).transformProperties
        (
            cpp.separation(patchFacei_)
        );
    }
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitProcessorPatch
(
    const processorPolyPatch& spp,
    TrackData& td
)
{}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitWallPatch
(
    const wallPolyPatch& spp,
    TrackData&
)
{}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitPatch
(
    const polyPatch& spp,
    TrackData&
)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleIO.C"

// ************************************************************************* //
