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
    Foam::solidParticle

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

SourceFiles
    solidParticleI.H
    solidParticle.C
    solidParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticle_H
#define solidParticle_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParticleCloud;

/*---------------------------------------------------------------------------*\
                           Class solidParticle Declaration
\*---------------------------------------------------------------------------*/

class solidParticle
:
    public Particle<solidParticle>
{
    // Private member data

        //- Diameter
        scalar d_;

        //- Velocity of parcel
        vector U_;


public:

    friend class Cloud<solidParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackData
    {
        //- Reference to the cloud containing this particle
        solidParticleCloud& spc_;

        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& nuInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;


    public:

        bool switchProcessor;
        bool keepParticle;


        // Constructors

            inline trackData
            (
                solidParticleCloud& spc,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& nuInterp,
                const vector& g
            );


        // Member functions

            inline solidParticleCloud& spc();

            inline const interpolationCellPoint<scalar>& rhoInterp() const;

            inline const interpolationCellPoint<vector>& UInterp() const;

            inline const interpolationCellPoint<scalar>& nuInterp() const;

            inline const vector& g() const;
    };


    // Constructors

        //- Construct from components
        inline solidParticle
        (
            const Cloud<solidParticle>& c,
            const vector& position,
            const label celli,
            const scalar m,
            const vector& U
        );

        //- Construct from Istream
        solidParticle
        (
            const Cloud<solidParticle>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<solidParticle> clone() const
        {
            return autoPtr<solidParticle>(new solidParticle(*this));
        }


    // Member Functions

        // Access

            //- Return diameter
            inline scalar d() const;

            //- Return velocity
            inline const vector& U() const;

            //- The nearest distance to a wall that
            //  the particle can be in the n direction
            inline scalar wallImpactDistance(const vector& n) const;


        // Tracking

            //- Move
            bool move(trackData&);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                solidParticle::trackData& td,
                const label patchI
            );

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions without trackData
            bool hitPatch
            (
                const polyPatch& p,
                int& td,
                const label patchI
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                solidParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch without trackData
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                solidParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            //- without trackData
            void hitWallPatch
            (
                const wallPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                solidParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a polyPatch
            //- without trackData
            void hitPatch
            (
                const polyPatch&,
                int&
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties
            (
                const tensor& T
            );

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            void transformProperties
            (
                const vector& separation
            );


    // I-O

        static void readFields(Cloud<solidParticle>& c);

        static void writeFields(const Cloud<solidParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const solidParticle&);
};


template<>
inline bool contiguous<solidParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
