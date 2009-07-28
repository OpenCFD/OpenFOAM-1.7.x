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
    Foam::searchableSurfaceWithGaps

Description
    searchableSurface using multiple slightly shifted underlying surfaces
    to make sure pierces don't go through gaps:
    - shift test vector with two small vectors (of size gap_) perpendicular
      to the original.
      Test with + and - this vector. Only if both register a hit is it seen
      as one.
    - extend the test vector slightly (with SMALL) to account for numerical
      inaccuracies.

SourceFiles
    searchableSurfaceWithGaps.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceWithGaps_H
#define searchableSurfaceWithGaps_H

#include "searchableSurface.H"
#include "UPtrList.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class searchableSurfaceWithGaps Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceWithGaps
:
    public searchableSurface
{
private:

    // Private Member Data

        //- gap size in meter
        const scalar gap_;

        //- Underlying geometry (size 1)
        UPtrList<searchableSurface> subGeom_;


    // Private Member Functions

        Pair<vector> offsetVecs(const point&, const point&) const;

        void offsetVecs
        (
            const pointField& start,
            const pointField& end,
            pointField& offset0,
            pointField& offset1
        ) const;

        static label countMisses
        (
            const List<pointIndexHit>& info,
            labelList& missMap
        );

        static label countMisses
        (
            const List<pointIndexHit>& plusInfo,
            const List<pointIndexHit>& minInfo,
            labelList& missMap
        );


        //- Disallow default bitwise copy construct
        searchableSurfaceWithGaps(const searchableSurfaceWithGaps&);

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaceWithGaps&);


public:

    //- Runtime type information
    TypeName("searchableSurfaceWithGaps");


    // Constructors

        //- Construct from dictionary (used by searchableSurface)
        searchableSurfaceWithGaps
        (
            const IOobject& io,
            const dictionary& dict
        );

    // Destructor

        virtual ~searchableSurfaceWithGaps();


    // Member Functions

        const searchableSurface& surface() const
        {
            return subGeom_[0];
        }


        virtual const wordList& regions() const
        {
            return surface().regions();
        }

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return surface().hasVolumeType();
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return surface().size();
        }


        // Multiple point queries.

            //- Find nearest on original surface. Note:does not use perturbation
            //  and hence might be inconsistent with intersections.
            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>& info
            ) const
            {
                surface().findNearest
                (
                    sample,
                    nearestDistSqr,
                    info
                );
            }

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit> >&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>& info,
                labelList& region
            ) const
            {
                surface().getRegion(info, region);
            }

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>& info,
                vectorField& normal
            ) const
            {
                surface().getNormal(info, normal);
            }

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField& samples,
                List<volumeType>& info
            ) const
            {
                surface().getVolumeType(samples, info);
            }


        // regIOobject implementation

            bool writeData(Ostream& os) const
            {
                return surface().writeData(os);
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
