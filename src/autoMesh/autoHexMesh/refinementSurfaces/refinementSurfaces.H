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
    Foam::refinementSurfaces

Description
    Container for data on surfaces used for surface-driven refinement.
    Contains all the data about the level of refinement needed per
    surface.

SourceFiles
    refinementSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef refinementSurfaces_H
#define refinementSurfaces_H

#include "triSurfaceGeoMesh.H"
#include "triSurfaceFields.H"
#include "vectorList.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;
class shellSurfaces;
class triSurfaceMesh;

/*---------------------------------------------------------------------------*\
                           Class refinementSurfaces Declaration
\*---------------------------------------------------------------------------*/

class refinementSurfaces
{
    // Private data

        //- Reference to all geometry.
        const searchableSurfaces& allGeometry_;

        //- Indices of surfaces that are refinement ones
        labelList surfaces_;

        //- Surface name (word)
        wordList names_;

        //- Per 'interface' surface : name of faceZone to put faces into
        wordList faceZoneNames_;

        //- Per 'interface' surface : name of cellZone to put cells into
        wordList cellZoneNames_;

        //- Per 'interface' surface : (only used if surface is closed)
        //  whether to zone cells inside or outside surface.
        boolList zoneInside_;


        //- From local region number to global region number
        labelList regionOffset_;

        //- From global region number to refinement level
        labelList minLevel_;

        //- From global region number to refinement level
        labelList maxLevel_;

        //- From global region number to perpendicular angle
        scalarField perpendicularAngle_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        refinementSurfaces(const refinementSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const refinementSurfaces&);


public:

    // Constructors

        //- Construct from surfaces and dictionaries
        refinementSurfaces
        (
            const searchableSurfaces& allGeometry,
            const PtrList<dictionary>&
        );

        //- Construct from surfaces and dictionary
        refinementSurfaces
        (
            const searchableSurfaces& allGeometry,
            const dictionary&
        );


    // Member Functions

        // Access

            const searchableSurfaces& geometry() const
            {
                return allGeometry_;
            }

            const labelList& surfaces() const
            {
                return surfaces_;
            }

            //- Names of surfaces
            const wordList& names() const
            {
                return names_;
            }

            //- Per 'interface' surface : name of faceZone to put faces into
            const wordList& faceZoneNames() const
            {
                return faceZoneNames_;
            }

            //- Per 'interface' surface : name of cellZone to put cells into
            const wordList& cellZoneNames() const
            {
                return cellZoneNames_;
            }

            //- Get indices of unnamed surfaces (surfaces without faceZoneName)
            labelList getUnnamedSurfaces() const;

            //- Get indices of named surfaces (surfaces with faceZoneName)
            labelList getNamedSurfaces() const;

            //- Get indices of closed named surfaces
            labelList getClosedNamedSurfaces() const;

            //- From local region number to global region number
            const labelList& regionOffset() const
            {
                return regionOffset_;
            }

            //- From global region number to refinement level
            const labelList& minLevel() const
            {
                return minLevel_;
            }

            //- From global region number to refinement level
            const labelList& maxLevel() const
            {
                return maxLevel_;
            }

            //- From global region number to perpendicular angle
            const scalarField& perpendicularAngle() const
            {
                return perpendicularAngle_;
            }


        // Helper

            //- From surface and region on surface to global region
            label globalRegion(const label surfI, const label regionI) const
            {
                return regionOffset_[surfI]+regionI;
            }

            //- Min level for surface and region on surface
            label minLevel(const label surfI, const label regionI) const
            {
                return minLevel_[globalRegion(surfI, regionI)];
            }

            //- Max level for surface and region on surface
            label maxLevel(const label surfI, const label regionI) const
            {
                return maxLevel_[globalRegion(surfI, regionI)];
            }

            label nRegions() const
            {
                return minLevel_.size();
            }

            //- Calculate minLevelFields
            void setMinLevelFields
            (
                const shellSurfaces& shells
            );

            //- Helper: count number of triangles per region
            static labelList countRegions(const triSurface&);


        // Searching

            //- Find intersection of edge. Return -1 or first surface
            //  with higher (than currentLevel) minlevel.
            //  Return surface number and level.
            void findHigherIntersection
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                labelList& surfaces,
                labelList& surfaceLevel
            ) const;

            //- Find all intersections of edge. Unsorted order.
            void findAllHigherIntersections
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                List<vectorList>& surfaceNormal,
                labelListList& surfaceLevel
            ) const;

            //- Find intersection nearest to the endpoints. surface1,2 are
            //  not indices into surfacesToTest but refinement surface indices.
            //  Returns surface, region on surface (so not global surface)
            //  and position on surface.
            void findNearestIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2
            ) const;

            //- Used for debugging only: find intersection of edge.
            void findAnyIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find nearest point on surfaces.
            void findNearest
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find nearest point on surfaces. Return surface and region on
            //  surface (so not global surface)
            void findNearestRegion
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const  scalarField& nearestDistSqr,
                labelList& hitSurface,
                labelList& hitRegion
            ) const;

            //- Detect if a point is 'inside' (closed) surfaces.
            //  Returns -1 if not, returns first surface it is.
            void findInside
            (
                const labelList& surfacesToTest,
                const pointField& pt,
                labelList& insideSurfaces
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
