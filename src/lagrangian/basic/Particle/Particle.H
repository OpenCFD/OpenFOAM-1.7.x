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
    Foam::Particle

Description

\*---------------------------------------------------------------------------*/

#ifndef Particle_H
#define Particle_H

#include "vector.H"
#include "IDLList.H"
#include "labelList.H"
#include "pointField.H"
#include "faceList.H"
#include "typeInfo.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Particle>
class Cloud;

class wedgePolyPatch;
class symmetryPolyPatch;
class cyclicPolyPatch;
class processorPolyPatch;
class wallPolyPatch;
class polyPatch;

// Forward declaration of friend functions and operators

template<class ParticleType>
class Particle;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const Particle<ParticleType>&
);


/*---------------------------------------------------------------------------*\
                           Class Particle Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class Particle
:
    public IDLList<ParticleType>::link
{

public:

    //- Class used to pass tracking data to the trackToFace function
    class trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            Cloud<ParticleType>& cloud_;


    public:

        bool switchProcessor;
        bool keepParticle;


        // Constructors

            inline trackData(Cloud<ParticleType>& cloud);


        // Member functions

            //- Return a reference to the cloud
            inline Cloud<ParticleType>& cloud();
    };


protected:

    // Private data

        //- Reference to the particle cloud
        const Cloud<ParticleType>& cloud_;

        //- Position of particle
        vector position_;

        //- Index of the cell it is in
        label celli_;

        //- Face index if the particle is on a face otherwise -1
        label facei_;

        //- Fraction of time-step completed
        scalar stepFraction_;

        //- Originating processor id
        label origProc_;

        //- Local particle id on originating processor
        label origId_;


    // Private member functions

        //- Return the 'lambda' value for the position, p, on the face,
        // where, p = from + lamda*(to - from)
        // for non-static meshes
        inline scalar lambda
        (
            const vector& from,
            const vector& to,
            const label facei,
            const scalar stepFraction
        ) const;

        //- Return the 'lambda' value for the position, p, on the face,
        // where, p = from + lamda*(to - from)
        // for static meshes
        inline scalar lambda
        (
            const vector& from,
            const vector& to,
            const label facei
        ) const;

        //- Find the faces between position and cell centre
        void findFaces
        (
            const vector& position,
            DynamicList<label>& faceList
        ) const;

        //- Find the faces between position and cell centre
        void findFaces
        (
            const vector& position,
            const label celli,
            const scalar stepFraction,
            DynamicList<label>& faceList
        ) const;


    // Patch interactions

        //- Overridable function to handle the particle hitting a patch
        //  Executed before other patch-hitting functions
        template<class TrackData>
        bool hitPatch
        (
            const polyPatch&,
            TrackData& td,
            const label patchI
        );

        //- Overridable function to handle the particle hitting a wedgePatch
        template<class TrackData>
        void hitWedgePatch
        (
            const wedgePolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a
        //  symmetryPatch
        template<class TrackData>
        void hitSymmetryPatch
        (
            const symmetryPolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a cyclicPatch
        template<class TrackData>
        void hitCyclicPatch
        (
            const cyclicPolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a
        //  processorPatch
        template<class TrackData>
        void hitProcessorPatch
        (
            const processorPolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        template<class TrackData>
        void hitWallPatch
        (
            const wallPolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a
        //  general patch
        template<class TrackData>
        void hitPatch
        (
            const polyPatch&,
            TrackData& td
        );


    // Transformations

        //- Transform the position the particle
        //  according to the given transformation tensor
        virtual void transformPosition(const tensor& T);

        //- Transform the physical properties of the particle
        //  according to the given transformation tensor
        virtual void transformProperties(const tensor& T);

        //- Transform the physical properties of the particle
        //  according to the given separation vector
        virtual void transformProperties(const vector& separation);


    // Parallel transfer

        //- Convert global addressing to the processor patch
        //  local equivalents
        template<class TrackData>
        void prepareForParallelTransfer(const label patchi, TrackData& td);

        //- Convert processor patch addressing to the global equivalents
        //  and set the celli to the face-neighbour
        template<class TrackData>
        void correctAfterParallelTransfer(const label patchi, TrackData& td);


public:

    friend class Cloud<ParticleType>;


    // Static data members

        //- Runtime type information
        TypeName("Particle");

        //- String representation of properties
        static string propHeader;


    // Constructors

        //- Construct from components
        Particle
        (
            const Cloud<ParticleType>&,
            const vector& position,
            const label celli
        );

        //- Construct from Istream
        Particle
        (
            const Cloud<ParticleType>&,
            Istream&,
            bool readFields = true
        );

        //- Construct as a copy
        Particle(const Particle& p);

        //- Construct a clone
        autoPtr<ParticleType> clone() const
        {
            return autoPtr<Particle>(new Particle(*this));
        }


        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            // Private data

                //- Reference to the cloud
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
                    new ParticleType
                    (
                        cloud_,
                        is,
                        true
                    )
                );
            }
        };


    //- Destructor
    virtual ~Particle()
    {}


    // Member Functions

        // Access

            //- Return true if particle is in cell
            inline bool inCell() const;

            //- Return true if position is in cell i
            inline bool inCell
            (
                const vector& position,
                const label celli,
                const scalar stepFraction
            ) const;

            //- Return current particle position
            inline const vector& position() const;

            //- Return current particle position
            inline vector& position();

            //- Return current cell particle is in
            inline label& cell();

            //- Return current cell particle is in
            inline label cell() const;

            //- Return current face particle is on otherwise -1
            inline label face() const;

            //- Return reference to the particle cloud
            inline const Cloud<ParticleType>& cloud() const;

            //- Return the impact model to be used, soft or hard (default).
            inline bool softImpact() const;

            //- Return the particle current time
            inline scalar currentTime() const;


        // Check

            //- Is the particle on the boundary/(or outside the domain)?
            inline bool onBoundary() const;

            //- Which patch is particle on
            inline label patch(const label facei) const;

            //- Which face of this patch is this particle on
            inline label patchFace
            (
                const label patchi,
                const label facei
            ) const;

            //- The nearest distance to a wall that
            //  the particle can be in the n direction
            inline scalar wallImpactDistance(const vector& n) const;

            //- Return the fraction of time-step completed
            inline scalar& stepFraction();

            //-  Return the fraction of time-step completed
            inline scalar stepFraction() const;

            //- Return the originating processor id
            inline label origProc() const;

            //- Return the particle id on originating processor
            inline label origId() const;


        // Track

            //- Track particle to end of trajectory
            //  or until it hits the boundary.
            //  On entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts and on exit it contains
            //  the fraction of the time-step completed.
            //  Returns the boundary face index if the track stops at the
            //  boundary, -1 otherwise.
            template<class TrackData>
            label track
            (
                const vector& endPosition,
                TrackData& td
            );

            //- Calls the templated track with dummy TrackData
            label track(const vector& endPosition);

            //- Track particle to a given position and returns 1.0 if the
            //  trajectory is completed without hitting a face otherwise
            //  stops at the face and returns the fraction of the trajectory
            //  completed.
            //  on entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts.
            template<class TrackData>
            scalar trackToFace
            (
                const vector& endPosition,
                TrackData& td
            );

            //- Calls the templated trackToFace with dummy TrackData
            scalar trackToFace(const vector& endPosition);

            //- Return the index of the face to be used in the interpolation
            //  routine
            inline label faceInterpolation() const;


    // I-O

        //- Read the fields associated with the owner cloud
        static void readFields(Cloud<ParticleType>& c);

        //- Write the fields associated with the owner cloud
        static void writeFields(const Cloud<ParticleType>& c);

        //- Write the particle data
        void write(Ostream& os, bool writeFields) const;

    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const Particle<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleI.H"

#define defineParticleTypeNameAndDebug(Type, DebugSwitch)                     \
    template<>                                                                \
    const Foam::word Particle<Type>::typeName(#Type);                         \
    template<>                                                                \
    int Particle<Type>::debug(Foam::debug::debugSwitch(#Type, DebugSwitch));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Particle.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
