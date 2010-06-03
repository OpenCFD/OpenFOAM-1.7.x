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
    Foam::syncTools

Description
    Various tools to aid synchronizing lists across coupled patches.

    Require
    - combineOperator (e.g. sumEqOp - not sumOp!) that is defined for the
      type and combineReduce(UList\<T\>, combineOperator) should be defined.
    - null value which gets overridden by any valid value.
    - transform function

    Can apply coordinate rotation/separation on cyclics but only for faces
    or if there is a single rotation/separation tensor.

SourceFiles
    syncTools.C
    syncToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef syncTools_H
#define syncTools_H

#include "UList.H"
#include "Pstream.H"
#include "transformList.H"
#include "Map.H"
#include "EdgeMap.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyBoundaryMesh;
class polyMesh;
class coupledPolyPatch;

/*---------------------------------------------------------------------------*\
                           Class syncTools Declaration
\*---------------------------------------------------------------------------*/

class syncTools
{
    // Private Member Functions

        //- Check whether uses couples.
        static bool hasCouples(const polyBoundaryMesh&);

        //- Check for single transformation tensor only.
        static void checkTransform(const coupledPolyPatch&, const bool);

        //- Apply separation to list. Either single vector or one vector
        //  per element.
        template <class T>
        static void separateList(const vectorField&, UList<T>&);

        template <class T>
        static void separateList(const vectorField&, Map<T>&);

        template <class T>
        static void separateList(const vectorField&, EdgeMap<T>&);

        //- Combine value with existing value in map.
        template <class T, class CombineOp>
        static void combine
        (
            Map<T>& pointValues,
            const CombineOp& cop,
            const label index,
            const T& val
        );

        //- Combine val with existing value at (implicit index) e.
        template <class T, class CombineOp>
        static void combine
        (
            EdgeMap<T>& edgeValues,
            const CombineOp& cop,
            const edge& index,
            const T& val
        );


public:

    // Static data members

        //- Synchronize values on all mesh points.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T, class CombineOp>
        static void syncPointList
        (
            const polyMesh&,
            UList<T>&,
            const CombineOp& cop,
            const T& nullValue,
            const bool applySeparation
        );

        //- Synchronize values on selected mesh points.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T, class CombineOp>
        static void syncPointList
        (
            const polyMesh&,
            const labelList& meshPoints,
            UList<T>&,
            const CombineOp& bop,
            const T& nullValue,
            const bool applySeparation
        );

        //- Synchronize values on all mesh edges.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T, class CombineOp>
        static void syncEdgeList
        (
            const polyMesh&,
            UList<T>&,
            const CombineOp& cop,
            const T& nullValue,
            const bool applySeparation
        );

        //- Synchronize values on boundary faces only.
        //  Optionally applies rotation tensor for non-parallel cyclics
        //  (but not separation!)
        template <class T, class CombineOp>
        static void syncBoundaryFaceList
        (
            const polyMesh&,
            UList<T>&,
            const CombineOp& cop,
            const bool applySeparation
        );

        //- Synchronize values on all mesh faces.
        //  Optionally applies rotation tensor for non-parallel cyclics
        //  (but not separation!)
        template <class T, class CombineOp>
        static void syncFaceList
        (
            const polyMesh&,
            UList<T>&,
            const CombineOp& cop,
            const bool applySeparation
        );

        //- Swap coupled face values.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T>
        static void swapBoundaryFaceList
        (
            const polyMesh&,
            UList<T>&,
            const bool applySeparation
        );

        //- Swap coupled face values.
        //  Applies rotation and optionally separation for parallel cyclics
        template <class T>
        static void swapFaceList
        (
            const polyMesh&,
            UList<T>&,
            const bool applySeparation
        );

        // Sparse versions

            //- Synchronize values on selected points.
            //  Applies rotation and optionally separation for parallel
            //  cyclics.
            template <class T, class CombineOp>
            static void syncPointMap
            (
                const polyMesh&,
                Map<T>& pointValues,
                const CombineOp& cop,
                const bool applySeparation
            );

            //- Synchronize values on selected edges. Edges are represented
            //  by the two vertices that make it up so global edges never get
            //  constructed.
            //  Applies rotation and optionally separation for parallel
            //  cyclics.
            template <class T, class CombineOp>
            static void syncEdgeMap
            (
                const polyMesh&,
                EdgeMap<T>& edgeValues,
                const CombineOp& cop,
                const bool applySeparation
            );

        // PackedList versions

            template <unsigned nBits, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                PackedList<nBits>& faceValues,
                const CombineOp& cop
            );

            template <unsigned nBits>
            static void swapFaceList
            (
                const polyMesh& mesh,
                PackedList<nBits>& faceValues
            );

            template <unsigned nBits, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                PackedList<nBits>& pointValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

            template <unsigned nBits, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                PackedList<nBits>& edgeValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

        // Other

            //- Get per point whether is it master (of a coupled set of points)
            static PackedBoolList getMasterPoints(const polyMesh&);

            //- Get per edge whether is it master (of a coupled set of edges)
            static PackedBoolList getMasterEdges(const polyMesh&);

            //- Get per face whether is it master (of a coupled set of faces)
            static PackedBoolList getMasterFaces(const polyMesh&);

};


template<>
void syncTools::separateList(const vectorField&, UList<vector>&);

template<>
void syncTools::separateList(const vectorField&, Map<vector>&);

template<>
void syncTools::separateList(const vectorField&, EdgeMap<vector>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "syncToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
