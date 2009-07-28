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
    Foam::searchableSurfaces

Description
    Container for searchableSurfaces.

SourceFiles
    searchableSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaces_H
#define searchableSurfaces_H

#include "searchableSurface.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class searchableSurfaces Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaces
:
    public PtrList<searchableSurface>
{
    // Private data

        //- Surface names
        wordList names_;

        //- Region names per surface
        List<wordList> regionNames_;

        ////- From global region name to surface and region on surface
        //HashTable<labelPair> regionNames_;

        //- Indices of all surfaces. Precalculated and stored.
        labelList allSurfaces_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        searchableSurfaces(const searchableSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaces&);


public:

    ClassName("searchableSurfaces");

    // Constructors

        //- Construct with length specified. Fill later.
        explicit searchableSurfaces(const label);


        ////- Construct from list of dictionaries
        //searchableSurfaces(const IOobject&, const PtrList<dictionary>&);

        //- Construct from dictionary
        searchableSurfaces(const IOobject&, const dictionary&);


    // Member Functions

        const wordList& names() const
        {
            return names_;
        }
        wordList& names()
        {
            return names_;
        }

        const List<wordList>& regionNames() const
        {
            return regionNames_;
        }
        List<wordList>& regionNames()
        {
            return regionNames_;
        }


        ////- If adding surfaces 'by hand'
        //HashTable<labelPair>& regionNames()
        //{
        //    return regionNames_;
        //}
        ////- Get surface and region for a name
        //const labelPair& surfaceRegion(const word& globalRegion) const
        //{
        //    return regionNames_[globalRegion];
        //}

        //- Find index of surface. Return -1 if not found.
        label findSurfaceID(const word& name) const;


        // Multiple point queries.

            //- Find any intersection. Return hit point information and
            //  surface number. If multiple surfaces hit the first surface
            //  is returned, not necessarily the nearest (to start).
            void findAnyIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find all intersections in order from start to end. Returns for
            //  every hit the surface and the hit info.
            void findAllIntersections
            (
                const pointField& start,
                const pointField& end,
                labelListList& surfaces,
                List<List<pointIndexHit> >&
            ) const;

            //- Find nearest. Return -1 (and a miss()) or surface and nearest
            //  point.
            void findNearest
            (
                const pointField&,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;


        // Single point queries

            //- Calculate point which is on a set of surfaces.
            pointIndexHit facesIntersection
            (
                const scalar initialDistSqr,
                const scalar convergenceDistSqr,
                const point& start
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
