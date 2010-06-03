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

#include "ExactParticle.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackingData>
Foam::label Foam::ExactParticle<ParticleType>::track
(
    const vector& endPosition,
    TrackingData& td
)
{
    this->facei_ = -1;

    // Tracks to endPosition or stop on boundary
    while (!this->onBoundary() && this->stepFraction_ < 1.0 - SMALL)
    {
        this->stepFraction_ +=
            trackToFace(endPosition, td)
           *(1.0 - this->stepFraction_);
    }

    return this->facei_;
}


template<class ParticleType>
Foam::label Foam::ExactParticle<ParticleType>::track
(
    const vector& endPosition
)
{
    int dummyTd;
    return track(endPosition, dummyTd);
}


template<class ParticleType>
template<class TrackingData>
Foam::scalar Foam::ExactParticle<ParticleType>::trackToFace
(
    const vector& endPosition,
    TrackingData& td
)
{
    const polyMesh& mesh = this->cloud().pMesh();
    const labelList& cFaces = mesh.cells()[this->celli_];

    point intersection(vector::zero);
    scalar trackFraction = VGREAT;
    label hitFacei = -1;

    const vector vec = endPosition-this->position_;

    forAll(cFaces, i)
    {
        label facei = cFaces[i];

        if (facei != this->face())
        {
            pointHit inter = mesh.faces()[facei].intersection
            (
                this->position_,
                vec,
                mesh.faceCentres()[facei],
                mesh.points(),
                intersection::HALF_RAY
            );

            if (inter.hit() && inter.distance() < trackFraction)
            {
                trackFraction = inter.distance();
                hitFacei = facei;
                intersection = inter.hitPoint();
            }
        }
    }

    if (hitFacei == -1)
    {
        // Did not find any intersection. Fall back to original approximate
        // algorithm
        return Particle<ParticleType>::trackToFace
        (
            endPosition,
            td
        );
    }

    if (trackFraction >= (1.0-SMALL))
    {
        // Nearest intersection beyond endPosition so we hit endPosition.
        trackFraction = 1.0;
        this->position_ = endPosition;
        this->facei_ = -1;
        return 1.0;
    }
    else
    {
        this->position_ = intersection;
        this->facei_ = hitFacei;
    }


    // Normal situation (trackFraction 0..1). Straight copy
    // of Particle::trackToFace.

    bool internalFace = this->cloud().internalFace(this->facei_);

    // change cell
    if (internalFace) // Internal face
    {
        if (this->celli_ == mesh.faceOwner()[this->facei_])
        {
            this->celli_ = mesh.faceNeighbour()[this->facei_];
        }
        else if (this->celli_ == mesh.faceNeighbour()[this->facei_])
        {
            this->celli_ = mesh.faceOwner()[this->facei_];
        }
        else
        {
            FatalErrorIn
            (
                "ExactParticle::trackToFace"
                "(const vector&, TrackingData&)"
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
            this->position_ = endPosition;
        }

        label patchi = patch(this->facei_);
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

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
        else if (isA<polyPatch>(patch))
        {
            p.hitPatch
            (
                static_cast<const polyPatch&>(patch), td
            );
        }
        else
        {
            FatalErrorIn
            (
                "ExactParticle::trackToFace"
                "(const vector& endPosition, scalar& trackFraction)"
            )<< "patch type " << patch.type() << " not suported" << nl
             << abort(FatalError);
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
        this->position_ +=
            1.0e-6*(mesh.cellCentres()[this->celli_] - this->position_);
    }

    return trackFraction;
}


template<class ParticleType>
Foam::scalar Foam::ExactParticle<ParticleType>::trackToFace
(
    const vector& endPosition
)
{
    int dummyTd;
    return trackToFace(endPosition, dummyTd);
}


template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ExactParticle<ParticleType>& p
)
{
    return operator<<(os, static_cast<const Particle<ParticleType>&>(p));
}


// ************************************************************************* //
