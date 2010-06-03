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
    Foam::trackedParticle

Description
    Particle class that marks cells it passes through. Used to mark cells
    visited by feature edges. Uses ExactParticle tracking class so
    will work on concave cells.

SourceFiles
    trackedParticle.C

\*---------------------------------------------------------------------------*/

#ifndef trackedParticle_H
#define trackedParticle_H

#include "ExactParticle.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class trackedParticleCloud;

/*---------------------------------------------------------------------------*\
                           Class trackedParticle Declaration
\*---------------------------------------------------------------------------*/

class trackedParticle
:
    public ExactParticle<trackedParticle>
{
    // Private data

        //- end point to track to
        point end_;

        //- level of this particle
        label level_;

        //- passive label
        label i_;

        //- passive label
        label j_;

public:

    friend class Cloud<trackedParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackData
    {
        //- Reference to the cloud containing this particle
        Cloud<trackedParticle>& cloud_;

        labelList& maxLevel_;

    public:

        bool switchProcessor;
        bool keepParticle;


        // Constructors

            trackData(Cloud<trackedParticle>& cloud, labelList& maxLevel)
            :
                cloud_(cloud),
                maxLevel_(maxLevel)
            {}


        // Member functions

            Cloud<trackedParticle>& cloud()
            {
                return cloud_;
            }

            labelList& maxLevel()
            {
                return maxLevel_;
            }
    };



    // Constructors

        //- Construct from components
        trackedParticle
        (
            const Cloud<trackedParticle>& c,
            const vector& position,
            const label celli,
            const point& end,
            const label level,
            const label i,
            const label j
        );

        //- Construct from Istream
        trackedParticle
        (
            const Cloud<trackedParticle>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<trackedParticle> clone() const
        {
            return autoPtr<trackedParticle>(new trackedParticle(*this));
        }


    // Member Functions

        //- point to track to
        point& end()
        {
            return end_;
        }

        //- transported label
        label& i()
        {
            return i_;
        }

        //- transported label
        label& j()
        {
            return j_;
        }



        // Tracking

            //- Track all particles to their end point
            bool move(trackData&);


            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                trackedParticle::trackData& td,
                const label patchI
            );
            bool hitPatch
            (
                const polyPatch&,
                int&,
                const label patchI
           );

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch
            (
                const wedgePolyPatch&,
                trackedParticle::trackData& td
            );
            void hitWedgePatch
            (
                const wedgePolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a
            //  symmetryPlane
            void hitSymmetryPatch
            (
                const symmetryPolyPatch&,
                trackedParticle::trackData& td
            );
            void hitSymmetryPatch
            (
                const symmetryPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch
            (
                const cyclicPolyPatch&,
                trackedParticle::trackData& td
            );
            void hitCyclicPatch
            (
                const cyclicPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackedParticle::trackData& td
            );
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                trackedParticle::trackData& td
            );
            void hitWallPatch
            (
                const wallPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackedParticle::trackData& td
            );
            void hitPatch
            (
                const polyPatch&,
                int&
            );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const trackedParticle&);


};


template<>
inline bool contiguous<trackedParticle>()
{
    return true;
}

//template<>
//void Cloud<trackedParticle>::readFields();
//
//template<>
//void Cloud<trackedParticle>::writeFields() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
