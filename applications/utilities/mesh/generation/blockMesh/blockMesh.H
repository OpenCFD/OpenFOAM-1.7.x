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
    Foam::blockMesh

Description

SourceFiles
    blockMesh.C

\*---------------------------------------------------------------------------*/

#ifndef blockMesh_H
#define blockMesh_H

#include "blockList.H"
#include "polyMesh.H"
#include "IOdictionary.H"
#include "curvedEdgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blockMesh Declaration
\*---------------------------------------------------------------------------*/

class blockMesh
:
    public blockList
{
    // Private data

        label nPoints_;
        label nCells_;

        curvedEdgeList edges_;

        polyMesh* topologyPtr_;

        labelList blockOffsets_;
        labelList mergeList_;

        pointField points_;
        cellShapeList cells_;
        faceListList patches_;


    // Private Member Functions

        bool blockLabelsOK
        (
            const label blockLabel,
            const pointField& points,
            const cellShape& blockShape
        );

        bool patchLabelsOK
        (
            const label patchLabel,
            const pointField& points,
            const faceList& patchShapes
        );

        polyMesh* createTopology(IOdictionary&);
        void checkBlockMesh(const polyMesh&);

        labelList createBlockOffsets();
        labelList createMergeList();

        pointField createPoints(const dictionary&);
        cellShapeList createCells();

        faceList createPatchFaces(const polyPatch& patchTopologyFaces);
        faceListList createPatches();

        //- as copy (not implemented)
        blockMesh(const blockMesh&);


public:

    // Constructors

        //- Construct from IOdictionary
        blockMesh(IOdictionary&);


    // Destructor

        ~blockMesh();


    // Member Functions

        // Access

            const polyMesh& topology() const;

            const curvedEdgeList& edges() const
            {
                return edges_;
            }

            const pointField& points() const
            {
                return points_;
            }

            const cellShapeList& cells() const
            {
                return cells_;
            }

            const faceListList& patches() const
            {
                return patches_;
            }

            wordList patchNames() const;

            wordList patchTypes() const;

            wordList patchPhysicalTypes() const;

            //- Number of blocks with specified zones
            label numZonedBlocks() const;


        // Write

            //- Writes edges of blockMesh in OBJ format.
            void writeTopology(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

