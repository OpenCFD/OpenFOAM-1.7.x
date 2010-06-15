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

Class
    Foam::ExactParticle

Description
    Special version of Particle to do tracking on non-convex cells.

\*---------------------------------------------------------------------------*/

#ifndef ExactParticle_H
#define ExactParticle_H

#include "face.H"
#include "Particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ExactParticle>
class Cloud;

// Forward declaration of friend functions and operators

template<class ParticleType>
class ExactParticle;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const ExactParticle<ParticleType>&
);


/*---------------------------------------------------------------------------*\
                           Class ExactParticle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class ExactParticle
:
    public Particle<ParticleType>
{

public:

    friend class Cloud<ParticleType>;


    // Constructors

        //- Construct from components
        ExactParticle
        (
            const Cloud<ParticleType>& cloud,
            const vector& position,
            const label celli
        )
        :
            Particle<ParticleType>(cloud, position, celli)
        {}


        //- Construct from Istream
        ExactParticle
        (
            const Cloud<ParticleType>& cloud,
            Istream& is,
            bool readFields = true
        )
        :
            Particle<ParticleType>(cloud, is, readFields)
        {}


        //- Factory class to read-construct particles used for parallel transfer
        class iNew
        {

            // Private data

            const Cloud<ParticleType>& cloud_;


        public:

            iNew(const Cloud<ParticleType>& cloud)
            :
                cloud_(cloud)
            {}

            autoPtr<ParticleType> operator()(Istream& is) const
            {
                return autoPtr<ParticleType>
                (
                    new ParticleType(cloud_, is)
                );
            }
        };


    // Destructor

        virtual ~ExactParticle()
        {}


    // Member Functions

            //- Track particle to end of trajectory
            //  or until it hits the boundary.
            //  On entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts and on exit it contains
            //  the fraction of the time-step completed.
            //  Returns the boundary face index if the track stops at the
            //  boundary, -1 otherwise.
            template<class TrackingData>
            label track
            (
                const vector& endPosition,
                TrackingData& td
            );

            //- Calls the templated track with dummy TrackingData
            label track(const vector& endPosition);

            //- Track particle to a given position and returns 1.0 if the
            //  trajectory is completed without hitting a face otherwise
            //  stops at the face and returns the fraction of the trajectory
            //  completed.
            //  on entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts.
            template<class TrackingData>
            scalar trackToFace
            (
                const vector& endPosition,
                TrackingData& td
            );

            //- Calls the templated trackToFace with dummy TrackingData
            scalar trackToFace(const vector& endPosition);


    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const ExactParticle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ExactParticle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
