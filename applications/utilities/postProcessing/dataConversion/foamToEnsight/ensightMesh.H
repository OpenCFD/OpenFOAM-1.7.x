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
    Foam::ensightMesh

Description

SourceFiles
    ensightMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ensightMesh_H
#define ensightMesh_H

#include "cellSets.H"
#include "faceSets.H"
#include "HashTable.H"
#include "HashSet.H"
#include "fvMesh.H"
#include "OFstream.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class argList;

/*---------------------------------------------------------------------------*\
                           Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightMesh
{
        class nFacePrimitives
        {
        public:

            label nPoints;
            label nTris;
            label nQuads;
            label nPolys;

            nFacePrimitives()
            :
                nPoints(0),
                nTris(0),
                nQuads(0),
                nPolys(0)
            {}
        };


    // Private data

        //- Reference to the OpenFOAM mesh
        const fvMesh& mesh_;

        //- Set binary file output
        bool binary_;

        //- The ensight part id for the first patch
        label patchPartOffset_;

        cellSets meshCellSets_;

        List<faceSets> boundaryFaceSets_;

        wordList allPatchNames_;

        List<labelList> allPatchProcs_;

        wordHashSet patchNames_;

        HashTable<nFacePrimitives> nPatchPrims_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ensightMesh(const ensightMesh&);

        //- Disallow default bitwise assignment
        void operator=(const ensightMesh&);

        void writePoints
        (
            const scalarField& pointsComponent,
            OFstream& ensightGeometryFile
        ) const;

        cellShapeList map
        (
            const cellShapeList& cellShapes,
            const labelList& prims
        ) const;

        cellShapeList map
        (
            const cellShapeList& cellShapes,
            const labelList& hexes,
            const labelList& wedges
        ) const;

        void writePrims
        (
            const cellShapeList& cellShapes,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        void writePolys
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        void writeAllPrims
        (
            const char* key,
            const label nPrims,
            const cellShapeList& cellShapes,
            const labelList& pointOffsets,
            OFstream& ensightGeometryFile
        ) const;

        void writeFacePrims
        (
            const char* key,
            const faceList& patchFaces,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        faceList map
        (
            const faceList& patchFaces,
            const labelList& prims
        ) const;

        void writeAllFacePrims
        (
            const char* key,
            const labelList& prims,
            const label nPrims,
            const faceList& patchFaces,
            const labelList& pointOffsets,
            const labelList& patchProcessors,
            OFstream& ensightGeometryFile
        ) const;

        void writeAscii
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;

        void writeBinary
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;

        void writePrimsBinary
        (
            const cellShapeList& cellShapes,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeAllPrimsBinary
        (
            const char* key,
            const label nPrims,
            const cellShapeList& cellShapes,
            const labelList& pointOffsets,
            std::ofstream& ensightGeometryFile
        ) const;

        void writePolysBinary
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeAllFacePrimsBinary
        (
            const char* key,
            const labelList& prims,
            const label nPrims,
            const faceList& patchFaces,
            const labelList& pointOffsets,
            const labelList& patchProcessors,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeFacePrimsBinary
        (
            const char* key,
            const faceList& patchFaces,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;


public:

    // Constructors

        //- Construct from fvMesh
        ensightMesh(const fvMesh&, const argList& args, const bool binary);


    //- Destructor
    ~ensightMesh();


    // Member Functions

        // Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }

            const cellSets& meshCellSets() const
            {
                return meshCellSets_;
            }

            const List<faceSets>& boundaryFaceSets() const
            {
                return boundaryFaceSets_;
            }

            const wordList& allPatchNames() const
            {
                return allPatchNames_;
            }

            const List<labelList>& allPatchProcs() const
            {
                return allPatchProcs_;
            }

            const wordHashSet& patchNames() const
            {
                return patchNames_;
            }

            const HashTable<nFacePrimitives>& nPatchPrims() const
            {
                return nPatchPrims_;
            }

            //- The ensight part id for the first patch
            label patchPartOffset() const
            {
                return patchPartOffset_;
            }


    // I-O

        void write
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
