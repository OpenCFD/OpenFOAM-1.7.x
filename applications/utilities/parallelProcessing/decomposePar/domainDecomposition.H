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
    Foam::domainDecomposition

Description
    Automatic domain decomposition class for FOAM meshes

SourceFiles
    domainDecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef domainDecomposition_H
#define domainDecomposition_H

#include "fvMesh.H"
#include "labelList.H"
#include "SLList.H"
#include "PtrList.H"
#include "point.H"
#include "Time.H"
#include "volFields.H"

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

/*---------------------------------------------------------------------------*\
                           Class domainDecomposition Declaration
\*---------------------------------------------------------------------------*/

class domainDecomposition
:
    public fvMesh
{
    // Private data

        //- Mesh decomposition control dictionary
        IOdictionary decompositionDict_;

        //- Number of processors in decomposition
        label nProcs_;

        //- Is the decomposition data to be distributed for each processor
        bool distributed_;

        //- Processor label for each cell
        labelList cellToProc_;

        //- Labels of points for each processor
        labelListList procPointAddressing_;

        //- Labels of faces for each processor
        // Note: Face turning index is stored as the sign on addressing
        // Only the processor boundary faces are affected: if the sign of the
        // index is negative, the processor face is the reverse of the
        // original face. In order to do this properly, all face
        // indices will be incremented by 1 and the decremented as
        // necessary t avoid the problem of face number zero having no
        // sign.
        labelListList procFaceAddressing_;

        //- Labels of cells for each processor
        labelListList procCellAddressing_;

        //- Original patch index for every processor patch
        labelListList procBoundaryAddressing_;

        //- Sizes for processor mesh patches
        // Excludes inter-processor boundaries
        labelListList procPatchSize_;

        //- Start indices for processor patches
        // Excludes inter-processor boundaries
        labelListList procPatchStartIndex_;

        //- Neighbour processor ID for inter-processor boundaries
        labelListList procNeighbourProcessors_;

        //- Sizes for inter-processor patches
        labelListList procProcessorPatchSize_;

        //- Start indices for inter-processor patches
        labelListList procProcessorPatchStartIndex_;

        //- List of globally shared point labels
        labelList globallySharedPoints_;

        //- Are there cyclic-parallel faces
        bool cyclicParallel_;


    // Private Member Functions

        void distributeCells();

        //- Mark all elements with value or -2 if occur twice
        static void mark
        (
            const labelList& zoneElems,
            const label zoneI,
            labelList& elementToZone
        );

public:

    // Constructors

        //- Construct from IOobject
        domainDecomposition(const IOobject& io);


    // Destructor

        ~domainDecomposition();


    // Member Functions

        //- Number of processor in decomposition
        label nProcs() const
        {
            return nProcs_;
        }

        //- Is the decomposition data to be distributed for each processor
        bool distributed() const
        {
            return distributed_;
        }

        //- Decompose mesh. Optionally remove zero-sized patches.
        void decomposeMesh(const bool filterEmptyPatches);

        //- Write decomposition
        bool writeDecomposition();

        //- Cell-processor decomposition labels
        const labelList& cellToProc() const
        {
            return cellToProc_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
