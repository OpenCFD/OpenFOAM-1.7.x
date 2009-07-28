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

Class
    Foam::directMappedPatchBase

Description
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

Note
    Storage is not optimal. It temporary collects all (patch)face centres
    on all processors to keep the addressing calculation simple.

SourceFiles
    directMappedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedPatchBase_H
#define directMappedPatchBase_H

#include "pointField.H"
#include "Tuple2.H"
#include "pointIndexHit.H"
#include "mapDistribute.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class directMappedPatchBase Declaration
\*---------------------------------------------------------------------------*/

class directMappedPatchBase
{

public:

        //- Mesh items to sample
        enum sampleMode
        {
            NEARESTCELL,
            NEARESTPATCHFACE,
            NEARESTFACE
        };

private:

    // Private data

        static const NamedEnum<sampleMode, 3> sampleModeNames_;

        //- Patch to sample
        const polyPatch& patch_;

        //- Region to sample
        const word sampleRegion_;

        //- What to sample
        const sampleMode mode_;

        //- Patch (only if NEARESTBOUNDARY)
        const word samplePatch_;

        //- Offset vector
        const vector offset_;

        //- Same region
        const bool sameRegion_;


        // Derived information

            //- Communication schedule:
            //  - Cells/faces to sample per processor
            //  - Patch faces to receive per processor
            //  - schedule
            mutable autoPtr<mapDistribute> mapPtr_;


    // Private Member Functions

        //- Collect single list of samples and originating processor+face.
        void collectSamples
        (
            pointField&,
            labelList& patchFaceProcs,
            labelList& patchFaces,
            pointField& patchFc
        ) const;

        //- Find cells/faces containing samples
        void findSamples
        (
            const pointField&,
            labelList& sampleProcs,     // processor containing sample
            labelList& sampleIndices,   // local index of cell/face
            pointField& sampleLocations // actual representative location
        ) const;

        //- Calculate matching
        void calcMapping() const;


    // Private class

        //- Private class for finding nearest
        //  - point+local index
        //  - sqr(distance)
        //  - processor
        typedef Tuple2<pointIndexHit, Tuple2<scalar, label> > nearInfo;

        class nearestEqOp
        {

        public:

            void operator()(nearInfo& x, const nearInfo& y) const
            {
                if (y.first().hit())
                {
                    if (!x.first().hit())
                    {
                        x = y;
                    }
                    else if (y.second().first() < x.second().first())
                    {
                        x = y;
                    }
                }
            }
        };


public:

    //- Runtime type information
    TypeName("directMappedPatchBase");


    // Constructors

        //- Construct from components
        directMappedPatchBase(const polyPatch&);

        //- Construct from dictionary
        directMappedPatchBase(const polyPatch&, const dictionary&);

        //- Construct as copy, resetting patch
        directMappedPatchBase(const polyPatch&, const directMappedPatchBase&);


    //- Destructor
    virtual ~directMappedPatchBase();


    // Member functions

        void clearOut();

        //- What to sample
        const sampleMode& mode() const
        {
            return mode_;
        }

        //- Region to sample
        const word& sampleRegion() const
        {
            return sampleRegion_;
        }

        //- Patch (only if NEARESTBOUNDARY)
        const word& samplePatch() const
        {
            return samplePatch_;
        }

        //- Offset vector (from patch faces to destination mesh objects)
        const vector& offset() const
        {
            return offset_;
        }

        //- Return reference to the parallel distribution map
        const mapDistribute& map() const
        {
            if (mapPtr_.empty())
            {
                calcMapping();
            }
            return mapPtr_();
        }

        //- Cached sampleRegion != mesh.name()
        bool sameRegion() const
        {
            return sameRegion_;
        }

        //- Get the region mesh
        const polyMesh& sampleMesh() const;

        //- Get the patch on the region
        const polyPatch& samplePolyPatch() const;

        //- Write as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
