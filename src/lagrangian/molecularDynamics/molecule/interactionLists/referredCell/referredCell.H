/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::referredCell

Description

SourceFiles
    referredCellI.H
    referredCell.C

\*---------------------------------------------------------------------------*/

#ifndef referredCell_H
#define referredCell_H

#include "vector.H"
#include "vectorList.H"
#include "tensor.H"
#include "transform.H"
#include "DynamicList.H"
#include "labelList.H"
#include "edgeList.H"
#include "polyMesh.H"

#include "referredMolecule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class referredCell Declaration
\*---------------------------------------------------------------------------*/

class referredCell
:
    public DynamicList<referredMolecule>
{
    // Private data

        label sourceProc_;

        label sourceCell_;

        //- Referred vertex positions
        vectorList vertexPositions_;

        edgeList edges_;

        labelListList faces_;

        vectorList faceCentres_;

        vectorList faceAreas_;

        labelList realCellsForInteraction_;

        vector offset_;

        tensor rotation_;


    // Private Member Functions

        void setConstructionData
        (
            const polyMesh& mesh,
            const label sourceCell
        );

        void locallyMapEdgeList
        (
            const labelList& points,
            const edgeList& sourceCellEdges
        );

        void locallyMapFaceList
        (
            const labelList& points,
            const labelListList& sourceCellFaces
        );

        vector referPosition(const vector& positionToRefer);

        vectorList referPositions(const vectorList& positionsToRefer);

        vector rotateVector(const vector& vectorToRotate);

        vectorList rotateVectors(const vectorList& vectorsToRotate);


public:

    // Constructors

        //- Construct null

        referredCell();

        //- Construct from components with external edge information
        referredCell
        (
            const polyMesh& mesh,
            const label sourceProc,
            const label sourceCell,
            const vector& offset,
            const tensor& rotation
        );

        //- Construct from components with existing local edge information
        referredCell
        (
            const label sourceProc,
            const label sourceCell,
            const vectorList& vertexPositions,
            const edgeList& localEdges,
            const labelListList& localFaces,
            const vectorList& faceCentres,
            const vectorList& faceAreas,
            const vector& offset,
            const tensor& rotation
        );

        //- Construct from pair of face centers (c) and plain
        //  face normals (n) (no need to make unit vectors or
        //  reverse one direction)
        //  Order of vectors important (S = source, D = Destination).
        //  External edge information.

        referredCell
        (
            const polyMesh& mesh,
            const label sourceProc,
            const label sourceCell,
            const vector& cS,
            const vector& cD,
            const vector& nS,
            const vector& nD
        );


    // Destructor

        virtual ~referredCell();


    // Member Functions

        //- Take this referredCell object that has already had it's transform
        //  calculated and refer it on again, retaining same source info.
        referredCell reRefer
        (
            const vector& cS,
            const vector& cD,
            const vector& nS,
            const vector& nD
        );

        //- Use internal transformatation values to transform the given
        //  postion to its new location.
        vector referPosition(const vector& positionToRefer) const;

        //- Use internal transformatation values to transform the given
        //  list of postions to their new locations.
        vectorList referPosition(const vectorList& positionsToRefer) const;

        //- Use internal transformatation values to rotate the given vector
        vector rotateVector(const vector& vectorToRotate) const;

        //- Use internal transformatation values to rotate the given
        // list of vectors
        vectorList rotateVectors(const vectorList& vectorsToRotate) const;

        //- referInMols clears, the stored list of referred mols and takes
        //  in a list of referred mols coming from a source processor,
        //  transformingtheir positions
        void referInMols(const List<referredMolecule>& incomingMols);

        //- duplicate() function to test whether a referred or real cell
        //  supplied by arguement is a duplicate of this referredCell.
        //  Can be used bi-directionally - i.e. can be called on an existing
        //  referred cell with a proposed referredCell as argument,
        //  or vice versa.  Can only be called by a proposed referredCell with
        //  a real cell index as arguement to test to see if the proposed
        //  referredCell is a duplicate.
        //  A duplicate cell is defined as one which has the same source
        //  processor, source cell, and an equal offset. Real cells have zero
        //  offset by definition.
        bool duplicate(const referredCell& refCellDupl) const;

        bool duplicate(const label procNo, const label nCells) const;


        // Access

            inline label sourceProc() const;

            inline label sourceCell() const;

            inline const vector& offset() const;

            inline const tensor& rotation() const;

            inline const vectorList& vertexPositions() const;

            inline const edgeList& edges() const;

            inline const labelListList& faces() const;

            inline const vectorList& faceCentres() const;

            inline const vectorList& faceAreas() const;

            inline labelList& realCells();

            inline const labelList& realCellsForInteraction() const;


    // Friend Operators

        inline friend bool operator==
        (
            const referredCell& a,
            const referredCell& b
        );

        inline friend bool operator!=
        (
            const referredCell& a,
            const referredCell& b
        );


    // IOstream Operators

        friend Istream& operator>>(Istream&, referredCell&);
        friend Ostream& operator<<(Ostream&, const referredCell&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "referredCellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
