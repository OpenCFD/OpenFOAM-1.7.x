/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::meshSearch

Description
    Various (local, not parallel) searches on polyMesh;
    uses (demand driven) octree to search.

SourceFiles
    meshSearch.C

\*---------------------------------------------------------------------------*/

#ifndef meshSearch_H
#define meshSearch_H

#include "pointIndexHit.H"
#include "Cloud.H"
#include "passiveParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class treeDataCell;
class treeDataFace;
class treeDataPoint;
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                           Class meshSearch Declaration
\*---------------------------------------------------------------------------*/

class meshSearch
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Whether to use face decomposition for all geometric tests
        const bool faceDecomp_;

        //- Dummy cloud to put particles on for tracking.
        Cloud<passiveParticle> cloud_;

        //- demand driven octrees

        mutable indexedOctree<treeDataFace>* boundaryTreePtr_;
        mutable indexedOctree<treeDataCell>* cellTreePtr_;
        mutable indexedOctree<treeDataPoint>* cellCentreTreePtr_;


    // Private Member Functions

        //- Updates nearestI, nearestDistSqr from any closer ones.
        static bool findNearer
        (
            const point& sample,
            const pointField& points,
            label& nearestI,
            scalar& nearestDistSqr
        );

        //- Updates nearestI, nearestDistSqr from any selected closer ones.
        static bool findNearer
        (
            const point& sample,
            const pointField& points,
            const labelList& indices,
            label& nearestI,
            scalar& nearestDistSqr
        );


        // Cells

            //- nearest cell centre using octree
            label findNearestCellTree(const point&) const;

            //- nearest cell centre going through all cells
            label findNearestCellLinear(const point&) const;

            //- walk from seed. Does not 'go around' boundary, just returns
            //  last cell before boundary.
            label findNearestCellWalk(const point&, const label) const;

            //- cell containing location. Linear search.
            label findCellLinear(const point&) const;


        // Cells

            label findNearestFaceTree(const point&) const;

            label findNearestFaceLinear(const point&) const;

            label findNearestFaceWalk(const point&, const label) const;



        // Boundary faces

            //- walk from seed to find nearest boundary face. Gets stuck in
            //  local minimum.
            label findNearestBoundaryFaceWalk
            (
                const point& location,
                const label seedFaceI
            ) const;

            //- Calculate offset vector in direction dir with as length a
            //  fraction of the cell size (of the cell straddling boundary face)
            vector offset
            (
                const point& bPoint,
                const label bFaceI,
                const vector& dir
            ) const;


        //- Disallow default bitwise copy construct
        meshSearch(const meshSearch&);

        //- Disallow default bitwise assignment
        void operator=(const meshSearch&);


public:

    // Declare name of the class and its debug switch
    ClassName("meshSearch");


    // Static data members

        //- tolerance on linear dimensions
        static scalar tol_;


    // Constructors

        //- Construct from components
        meshSearch(const polyMesh& mesh, const bool faceDecomp = true);


    // Destructor

        ~meshSearch();


    // Member Functions

        // Access

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            //- Get (demand driven) reference to octree holding all
            //  boundary faces
            const indexedOctree<treeDataFace>& boundaryTree() const;

            //- Get (demand driven) reference to octree holding all cells
            const indexedOctree<treeDataCell>& cellTree() const;

            //- Get (demand driven) reference to octree holding all cell centres
            const indexedOctree<treeDataPoint>& cellCentreTree() const;


        // Queries

            //- test for point in cell. Does not handle cells with center
            //  outside cell.
            bool pointInCell(const point& p, const label celli) const;

            //- Find nearest cell in terms of cell centre.
            // - use octree
            // - use linear search
            // - if seed is provided walk. (uses findNearestCellWalk; 
            //   does not handle holes in domain)
            label findNearestCell
            (
                const point& location,
                const label seedCellI = -1,
                const bool useTreeSearch = true
            ) const;

            label findNearestFace
            (
                const point& location,
                const label seedFaceI = -1,
                const bool useTreeSearch = true
            ) const;

            //- Find cell containing (using pointInCell) location.
            //  If seed provided walks and falls back to linear/tree search.
            //  (so handles holes correctly)s
            //  Returns -1 if not in domain.
            label findCell
            (
                const point& location,
                const label seedCellI = -1,
                const bool useTreeSearch = true
            ) const;

            //- Find nearest boundary face
            //  If seed provided walks but then does not pass local minima
            //  in distance. Also does not jump from one connected region to
            //  the next.
            label findNearestBoundaryFace
            (
                const point& location,
                const label seedFaceI = -1,
                const bool useTreeSearch = true
            ) const;

            //- Find first intersection of boundary in segment [pStart, pEnd]
            //  (so inclusive of endpoints). Always octree for now
            pointIndexHit intersection(const point& pStart, const point& pEnd)
            const;

            //- Find all intersections of boundary within segment pStart .. pEnd
            //  Always octree for now
            List<pointIndexHit> intersections
            (
                const point& pStart,
                const point& pEnd
            ) const;

            //- Determine inside/outside status
            bool isInside(const point&) const;


        //- delete all storage
        void clearOut();

        //- Correct for mesh geom/topo changes
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
