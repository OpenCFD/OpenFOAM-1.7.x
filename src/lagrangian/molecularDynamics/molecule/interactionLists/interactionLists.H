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
    Foam::interactionLists

Description

SourceFiles
    interactionListsI.H
    interactionLists.C
    interactionListsIO.C

\*---------------------------------------------------------------------------*/

#ifndef interactionLists_H
#define interactionLists_H

#include "polyMesh.H"
#include "vector2D.H"
#include "directInteractionList.H"
#include "referredCell.H"
#include "referredCellList.H"
#include "sendingReferralList.H"
#include "receivingReferralList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interactionLists Declaration
\*---------------------------------------------------------------------------*/

class interactionLists
{
    // Private data

        const polyMesh& mesh_;

        scalar rCutMaxSqr_;

        directInteractionList dil_;

        referredCellList ril_;

        List<sendingReferralList> cellSendingReferralLists_;

        List<receivingReferralList> cellReceivingReferralLists_;


    // Private Member Functions

        //- Build referralLists which define how to send information
        // to referredCells to source cells
        void buildCellReferralLists();

        //- Disallow default bitwise copy construct
        interactionLists(const interactionLists&);

        //- Disallow default bitwise assignment
        void operator=(const interactionLists&);


public:

    // Static data members

        //- Tolerance for checking that faces on a patch segment
        static scalar transTol;


    // Constructors

        //- Construct and create all information from the mesh
        interactionLists
        (
            const polyMesh& mesh,
            scalar rCutMaxSqr,
            bool pointPointListBuild = false
        );

        //- Construct from file
        interactionLists(const polyMesh& mesh);


    // Destructor

        ~interactionLists();


    // Member Functions

            bool testPointPointDistance
            (
                const label ptI,
                const label ptJ
            ) const;

            bool testPointFaceDistance
            (
                const label p,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const label p,
                const referredCell& refCell
            ) const;

            bool testPointFaceDistance
            (
                const vectorList& pointsToTest,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const vector& p,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const vector& p,
                const labelList& faceToTest,
                const vectorList& points,
                const vector& faceC,
                const vector& faceA
            ) const;

            bool testEdgeEdgeDistance
            (
                const edge& eI,
                const edge& eJ
            ) const;

            bool testEdgeEdgeDistance
            (
                const edge& eI,
                const vector& eJs,
                const vector& eJe
            ) const;

            const labelList realCellsInRangeOfSegment
            (
                const labelList& segmentFaces,
                const labelList& segmentEdges,
                const labelList& segmentPoints
            ) const;

            const labelList referredCellsInRangeOfSegment
            (
                const List<referredCell>& referredInteractionList,
                const labelList& segmentFaces,
                const labelList& segmentEdges,
                const labelList& segmentPoints
            ) const;


        // Access

            inline const polyMesh& mesh() const;

            inline const directInteractionList& dil() const;

            inline const referredCellList& ril() const;
            inline referredCellList& ril();

            inline const List<sendingReferralList>&
                cellSendingReferralLists() const;

            inline const List<receivingReferralList>&
                cellReceivingReferralLists() const;

            inline label nInteractingProcs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interactionListsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
