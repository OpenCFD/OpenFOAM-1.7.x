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
    Foam::searchableSurfaceCollection

Description
    Set of transformed searchableSurfaces. Does not do boolean operations.
    So when meshing might find parts 'inside'.

SourceFiles
    searchableSurfaceCollection.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceCollection_H
#define searchableSurfaceCollection_H

#include "searchableSurface.H"
#include "treeBoundBox.H"
#include "coordinateSystem.H"
#include "UPtrList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                 Class searchableSurfaceCollection Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceCollection
:
    public searchableSurface
{
private:

    // Private Member Data

        // Per instance data

            //- instance name
            wordList instance_;

            //- scaling vector
            vectorField scale_;

            //- transformation
            PtrList<coordinateSystem> transform_;

            UPtrList<searchableSurface> subGeom_;

            Switch mergeSubRegions_;

        //- Region names
        mutable wordList regions_;
        //- From individual regions to collection regions
        mutable labelList regionOffset_;


    // Private Member Functions

        //- Find point nearest to sample. Updates minDistSqr. Sets nearestInfo
        //  and surface index
        void findNearest
        (
            const pointField& samples,
            scalarField& minDistSqr,
            List<pointIndexHit>& nearestInfo,
            labelList& nearestSurf
        ) const;


        //- Disallow default bitwise copy construct
        searchableSurfaceCollection(const searchableSurfaceCollection&);

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaceCollection&);


public:

    //- Runtime type information
    TypeName("searchableSurfaceCollection");


    // Constructors

        //- Construct from dictionary (used by searchableSurface)
        searchableSurfaceCollection
        (
            const IOobject& io,
            const dictionary& dict
        );

    // Destructor

        virtual ~searchableSurfaceCollection();


    // Member Functions

        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return false;
        }

        //- Range of local indices that can be returned.
        virtual label size() const;


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

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
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                notImplemented
                (
                    "searchableSurfaceCollection::writeData(Ostream&) const"
                );
                return false;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
