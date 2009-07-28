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

\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "parMetisDecomp.H"
#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "polyMesh.H"
#include "Time.H"
#include "labelIOField.H"
#include "globalIndex.H"

#include <mpi.h>

extern "C"
{
#   include "parmetis.h"
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parMetisDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        parMetisDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Does prevention of 0 cell domains and calls parmetis.
Foam::label Foam::parMetisDecomp::decompose
(
    Field<int>& xadj,
    Field<int>& adjncy,
    const pointField& cellCentres,
    Field<int>& cellWeights,
    Field<int>& faceWeights,
    const List<int>& options,

    List<int>& finalDecomp
)
{
    // C style numbering
    int numFlag = 0;

    // Number of dimensions
    int nDims = 3;


    if (cellCentres.size() != xadj.size()-1)
    {
        FatalErrorIn("parMetisDecomp::decompose(..)")
            << "cellCentres:" << cellCentres.size()
            << " xadj:" << xadj.size()
            << abort(FatalError);
    }


    // Get number of cells on all processors
    List<int> nLocalCells(Pstream::nProcs());
    nLocalCells[Pstream::myProcNo()] = xadj.size()-1;
    Pstream::gatherList(nLocalCells);
    Pstream::scatterList(nLocalCells);

    // Get cell offsets.
    List<int> cellOffsets(Pstream::nProcs()+1);
    int nGlobalCells = 0;
    forAll(nLocalCells, procI)
    {
        cellOffsets[procI] = nGlobalCells;
        nGlobalCells += nLocalCells[procI];
    }
    cellOffsets[Pstream::nProcs()] = nGlobalCells;

    // Convert pointField into float
    Field<floatScalar> xyz(3*cellCentres.size());
    int compI = 0;
    forAll(cellCentres, cellI)
    {
        const point& cc = cellCentres[cellI];
        xyz[compI++] = float(cc.x());
        xyz[compI++] = float(cc.y());
        xyz[compI++] = float(cc.z());
    }

    // Make sure every domain has at least one cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (Metis falls over with zero sized domains)
    // Trickle cells from processors that have them up to those that
    // don't.


    // Number of cells to send to the next processor
    // (is same as number of cells next processor has to receive)
    List<int> nSendCells(Pstream::nProcs(), 0);

    for (label procI = nLocalCells.size()-1; procI >=1; procI--)
    {
        if (nLocalCells[procI]-nSendCells[procI] < 1)
        {
            nSendCells[procI-1] = nSendCells[procI]-nLocalCells[procI]+1;
        }
    }

    // First receive (so increasing the sizes of all arrays)

    if (Pstream::myProcNo() >= 1 && nSendCells[Pstream::myProcNo()-1] > 0)
    {
        // Receive cells from previous processor
        IPstream fromPrevProc(Pstream::blocking, Pstream::myProcNo()-1);

        Field<int> prevXadj(fromPrevProc);
        Field<int> prevAdjncy(fromPrevProc);
        Field<floatScalar> prevXyz(fromPrevProc);
        Field<int> prevCellWeights(fromPrevProc);
        Field<int> prevFaceWeights(fromPrevProc);

        if (prevXadj.size() != nSendCells[Pstream::myProcNo()-1])
        {
            FatalErrorIn("parMetisDecomp::decompose(..)")
                << "Expected from processor " << Pstream::myProcNo()-1
                << " connectivity for " << nSendCells[Pstream::myProcNo()-1]
                << " nCells but only received " << prevXadj.size()
                << abort(FatalError);
        }

        // Insert adjncy
        prepend(prevAdjncy, adjncy);
        // Adapt offsets and prepend xadj
        xadj += prevAdjncy.size();
        prepend(prevXadj, xadj);
        // Coords
        prepend(prevXyz, xyz);
        // Weights
        prepend(prevCellWeights, cellWeights);
        prepend(prevFaceWeights, faceWeights);
    }


    // Send to my next processor

    if (nSendCells[Pstream::myProcNo()] > 0)
    {
        // Send cells to next processor
        OPstream toNextProc(Pstream::blocking, Pstream::myProcNo()+1);

        int nCells = nSendCells[Pstream::myProcNo()];
        int startCell = xadj.size()-1 - nCells;
        int startFace = xadj[startCell];
        int nFaces = adjncy.size()-startFace;

        // Send for all cell data: last nCells elements
        // Send for all face data: last nFaces elements
        toNextProc
            << Field<int>::subField(xadj, nCells, startCell)-startFace
            << Field<int>::subField(adjncy, nFaces, startFace)
            << SubField<floatScalar>(xyz, nDims*nCells, nDims*startCell)
            <<
            (
                cellWeights.size()
              ? static_cast<const Field<int>&>
                (
                    Field<int>::subField(cellWeights, nCells, startCell)
                )
              : Field<int>(0)
            )
            <<
            (
                faceWeights.size()
              ? static_cast<const Field<int>&>
                (
                    Field<int>::subField(faceWeights, nFaces, startFace)
                )
              : Field<int>(0)
            );

        // Remove data that has been sent
        if (faceWeights.size())
        {
            faceWeights.setSize(faceWeights.size()-nFaces);
        }
        if (cellWeights.size())
        {
            cellWeights.setSize(cellWeights.size()-nCells);
        }
        xyz.setSize(xyz.size()-nDims*nCells);
        adjncy.setSize(adjncy.size()-nFaces);
        xadj.setSize(xadj.size() - nCells);
    }



    // Adapt number of cells
    forAll(nSendCells, procI)
    {
        // Sent cells
        nLocalCells[procI] -= nSendCells[procI];

        if (procI >= 1)
        {
            // Received cells
            nLocalCells[procI] += nSendCells[procI-1];
        }
    }
    // Adapt cellOffsets
    nGlobalCells = 0;
    forAll(nLocalCells, procI)
    {
        cellOffsets[procI] = nGlobalCells;
        nGlobalCells += nLocalCells[procI];
    }


    if (nLocalCells[Pstream::myProcNo()] != (xadj.size()-1))
    {
        FatalErrorIn("parMetisDecomp::decompose(..)")
            << "Have connectivity for " << xadj.size()-1
            << " cells but nLocalCells:" << nLocalCells[Pstream::myProcNo()]
            << abort(FatalError);
    }

    // Weight info
    int wgtFlag = 0;
    int* vwgtPtr = NULL;
    int* adjwgtPtr = NULL;

    if (cellWeights.size())
    {
        vwgtPtr = cellWeights.begin();
        wgtFlag += 2;       // Weights on vertices
    }
    if (faceWeights.size())
    {
        adjwgtPtr = faceWeights.begin();
        wgtFlag += 1;       // Weights on edges
    }


    // Number of weights or balance constraints
    int nCon = 1;
    // Per processor, per constraint the weight
    Field<floatScalar> tpwgts(nCon*nProcessors_, 1./nProcessors_);
    // Imbalance tolerance
    Field<floatScalar> ubvec(nCon, 1.02);
    if (nProcessors_ == 1)
    {
        // If only one processor there is no imbalance.
        ubvec[0] = 1;
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    // output: cell -> processor addressing
    finalDecomp.setSize(nLocalCells[Pstream::myProcNo()]);

    // output: number of cut edges
    int edgeCut = 0;


    ParMETIS_V3_PartGeomKway
    (
        cellOffsets.begin(),    // vtxDist
        xadj.begin(),
        adjncy.begin(),
        vwgtPtr,                // vertexweights
        adjwgtPtr,              // edgeweights
        &wgtFlag,
        &numFlag,
        &nDims,
        xyz.begin(),
        &nCon,
        &nProcessors_,          // nParts
        tpwgts.begin(),
        ubvec.begin(),
        const_cast<List<int>&>(options).begin(),
        &edgeCut,
        finalDecomp.begin(),
        &comm
    );


    // If we sent cells across make sure we undo it
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Receive back from next processor if I sent something
    if (nSendCells[Pstream::myProcNo()] > 0)
    {
        IPstream fromNextProc(Pstream::blocking, Pstream::myProcNo()+1);

        List<int> nextFinalDecomp(fromNextProc);

        if (nextFinalDecomp.size() != nSendCells[Pstream::myProcNo()])
        {
            FatalErrorIn("parMetisDecomp::decompose(..)")
                << "Expected from processor " << Pstream::myProcNo()+1
                << " decomposition for " << nSendCells[Pstream::myProcNo()]
                << " nCells but only received " << nextFinalDecomp.size()
                << abort(FatalError);
        }

        append(nextFinalDecomp, finalDecomp);
    }

    // Send back to previous processor.
    if (Pstream::myProcNo() >= 1 && nSendCells[Pstream::myProcNo()-1] > 0)
    {
        OPstream toPrevProc(Pstream::blocking, Pstream::myProcNo()-1);

        int nToPrevious = nSendCells[Pstream::myProcNo()-1];

        toPrevProc <<
            SubList<int>
            (
                finalDecomp,
                nToPrevious,
                finalDecomp.size()-nToPrevious
            );

        // Remove locally what has been sent
        finalDecomp.setSize(finalDecomp.size()-nToPrevious);
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parMetisDecomp::parMetisDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::parMetisDecomp::decompose(const pointField& points)
{
    if (points.size() != mesh_.nCells())
    {
        FatalErrorIn("parMetisDecomp::decompose(const pointField&)")
            << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    // For running sequential ...
    if (Pstream::nProcs() <= 1)
    {
        return metisDecomp(decompositionDict_, mesh_).decompose(points);
    }

    // Create global cell numbers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Get number of cells on all processors
    List<int> nLocalCells(Pstream::nProcs());
    nLocalCells[Pstream::myProcNo()] = mesh_.nCells();
    Pstream::gatherList(nLocalCells);
    Pstream::scatterList(nLocalCells);

    // Get cell offsets.
    List<int> cellOffsets(Pstream::nProcs()+1);
    int nGlobalCells = 0;
    forAll(nLocalCells, procI)
    {
        cellOffsets[procI] = nGlobalCells;
        nGlobalCells += nLocalCells[procI];
    }
    cellOffsets[Pstream::nProcs()] = nGlobalCells;

    int myOffset = cellOffsets[Pstream::myProcNo()];


    //
    // Make Metis Distributed CSR (Compressed Storage Format) storage
    //   adjncy      : contains cellCells (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    //



    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    // Get renumbered owner on other side of coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<int> globalNeighbour(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start() - mesh_.nInternalFaces();

            forAll(pp, i)
            {
                globalNeighbour[bFaceI++] = faceOwner[faceI++] + myOffset;
            }
        }
    }

    // Get the cell on the other side of coupled patches
    syncTools::swapBoundaryFaceList(mesh_, globalNeighbour, false);


    // Count number of faces (internal + coupled)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Number of faces per cell
    List<int> nFacesPerCell(mesh_.nCells(), 0);

    // Number of coupled faces
    label nCoupledFaces = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        nFacesPerCell[faceOwner[faceI]]++;
        nFacesPerCell[faceNeighbour[faceI]]++;
    }
    // Handle coupled faces
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                nCoupledFaces++;
                nFacesPerCell[faceOwner[faceI++]]++;
            }
        }
    }


    // Fill in xadj
    // ~~~~~~~~~~~~

    Field<int> xadj(mesh_.nCells()+1, -1);

    int freeAdj = 0;

    for (label cellI = 0; cellI < mesh_.nCells(); cellI++)
    {
        xadj[cellI] = freeAdj;

        freeAdj += nFacesPerCell[cellI];
    }
    xadj[mesh_.nCells()] = freeAdj;



    // Fill in adjncy
    // ~~~~~~~~~~~~~~

    Field<int> adjncy(2*mesh_.nInternalFaces() + nCoupledFaces, -1);

    nFacesPerCell = 0;

    // For internal faces is just offsetted owner and neighbour
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = faceOwner[faceI];
        label nei = faceNeighbour[faceI];

        adjncy[xadj[own] + nFacesPerCell[own]++] = nei + myOffset;
        adjncy[xadj[nei] + nFacesPerCell[nei]++] = own + myOffset;
    }
    // For boundary faces is offsetted coupled neighbour
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start()-mesh_.nInternalFaces();

            forAll(pp, i)
            {
                label own = faceOwner[faceI];
                adjncy[xadj[own] + nFacesPerCell[own]++] =
                    globalNeighbour[bFaceI];

                faceI++;
                bFaceI++;
            }
        }
    }



    // decomposition options. 0 = use defaults
    List<int> options(3, 0);
    //options[0] = 1;     // don't use defaults but use values below
    //options[1] = -1;    // full debug info
    //options[2] = 15;    // random number seed

    // cell weights (so on the vertices of the dual)
    Field<int> cellWeights;

    // face weights (so on the edges of the dual)
    Field<int> faceWeights;

    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("metisCoeffs"))
    {
        const dictionary& metisCoeffs =
            decompositionDict_.subDict("metisCoeffs");
        word weightsFile;

        if (metisCoeffs.readIfPresent("cellWeightsFile", weightsFile))
        {
            Info<< "parMetisDecomp : Using cell-based weights read from "
                << weightsFile << endl;

            labelIOField cellIOWeights
            (
                IOobject
                (
                    weightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
            cellWeights.transfer(cellIOWeights);

            if (cellWeights.size() != mesh_.nCells())
            {
                FatalErrorIn("parMetisDecomp::decompose(const pointField&)")
                    << "Number of cell weights " << cellWeights.size()
                    << " read from " << cellIOWeights.objectPath()
                    << " does not equal number of cells " << mesh_.nCells()
                    << exit(FatalError);
            }
        }

        if (metisCoeffs.readIfPresent("faceWeightsFile", weightsFile))
        {
            Info<< "parMetisDecomp : Using face-based weights read from "
                << weightsFile << endl;

            labelIOField weights
            (
                IOobject
                (
                    weightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );

            if (weights.size() != mesh_.nFaces())
            {
                FatalErrorIn("parMetisDecomp::decompose(const pointField&)")
                    << "Number of face weights " << weights.size()
                    << " does not equal number of internal and boundary faces "
                    << mesh_.nFaces()
                    << exit(FatalError);
            }

            faceWeights.setSize(2*mesh_.nInternalFaces()+nCoupledFaces);

            // Assume symmetric weights. Keep same ordering as adjncy.
            nFacesPerCell = 0;

            // Handle internal faces
            for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
            {
                label w = weights[faceI];

                label own = faceOwner[faceI];
                label nei = faceNeighbour[faceI];

                faceWeights[xadj[own] + nFacesPerCell[own]++] = w;
                faceWeights[xadj[nei] + nFacesPerCell[nei]++] = w;
            }
            // Coupled boundary faces
            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (pp.coupled())
                {
                    label faceI = pp.start();

                    forAll(pp, i)
                    {
                        label w = weights[faceI];
                        label own = faceOwner[faceI];
                        adjncy[xadj[own] + nFacesPerCell[own]++] = w;
                        faceI++;
                    }
                }
            }
        }

        if (metisCoeffs.readIfPresent("options", options))
        {
            Info<< "Using Metis options     " << options
                << nl << endl;

            if (options.size() != 3)
            {
                FatalErrorIn("parMetisDecomp::decompose(const pointField&)")
                    << "Number of options " << options.size()
                    << " should be three." << exit(FatalError);
            }
        }
    }


    // Do actual decomposition
    List<int> finalDecomp;
    decompose
    (
        xadj,
        adjncy,
        points,
        cellWeights,
        faceWeights,
        options,

        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


Foam::labelList Foam::parMetisDecomp::decompose
(
    const labelList& cellToRegion,
    const pointField& regionPoints
)
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (cellToRegion.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "parMetisDecomp::decompose(const labelList&, const pointField&)"
        )   << "Size of cell-to-coarse map " << cellToRegion.size()
            << " differs from number of cells in mesh " << mesh_.nCells()
            << exit(FatalError);
    }


    // Global region numbering engine
    globalIndex globalRegions(regionPoints.size());


    // Get renumbered owner region on other side of coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<int> globalNeighbour(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start() - mesh_.nInternalFaces();

            forAll(pp, i)
            {
                label ownRegion = cellToRegion[faceOwner[faceI]];
                globalNeighbour[bFaceI++] = globalRegions.toGlobal(ownRegion);
                faceI++;
            }
        }
    }

    // Get the cell on the other side of coupled patches
    syncTools::swapBoundaryFaceList(mesh_, globalNeighbour, false);


    // Get globalCellCells on coarse mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList globalRegionRegions;
    {
        List<DynamicList<label> > dynRegionRegions(regionPoints.size());

        // Internal faces first
        forAll(faceNeighbour, faceI)
        {
            label ownRegion = cellToRegion[faceOwner[faceI]];
            label neiRegion = cellToRegion[faceNeighbour[faceI]];

            if (ownRegion != neiRegion)
            {
                label globalOwn = globalRegions.toGlobal(ownRegion);
                label globalNei = globalRegions.toGlobal(neiRegion);

                if (findIndex(dynRegionRegions[ownRegion], globalNei) == -1)
                {
                    dynRegionRegions[ownRegion].append(globalNei);
                }
                if (findIndex(dynRegionRegions[neiRegion], globalOwn) == -1)
                {
                    dynRegionRegions[neiRegion].append(globalOwn);
                }
            }
        }

        // Coupled boundary faces
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                label faceI = pp.start();
                label bFaceI = pp.start() - mesh_.nInternalFaces();

                forAll(pp, i)
                {
                    label ownRegion = cellToRegion[faceOwner[faceI]];
                    label globalNei = globalNeighbour[bFaceI++];
                    faceI++;

                    if (findIndex(dynRegionRegions[ownRegion], globalNei) == -1)
                    {
                        dynRegionRegions[ownRegion].append(globalNei);
                    }
                }
            }
        }

        globalRegionRegions.setSize(dynRegionRegions.size());
        forAll(dynRegionRegions, i)
        {
            globalRegionRegions[i].transfer(dynRegionRegions[i]);
        }
    }

    labelList regionDecomp(decompose(globalRegionRegions, regionPoints));

    // Rework back into decomposition for original mesh_
    labelList cellDistribution(cellToRegion.size());

    forAll(cellDistribution, cellI)
    {
        cellDistribution[cellI] = regionDecomp[cellToRegion[cellI]];
    }
    return cellDistribution;
}


Foam::labelList Foam::parMetisDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "parMetisDecomp::decompose(const labelListList&, const pointField&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }

    // For running sequential ...
    if (Pstream::nProcs() <= 1)
    {
        return metisDecomp(decompositionDict_, mesh_)
            .decompose(globalCellCells, cellCentres);
    }


    // Make Metis Distributed CSR (Compressed Storage Format) storage

    // Connections
    Field<int> adjncy;
    // Offsets into adjncy
    Field<int> xadj;
    metisDecomp::calcMetisCSR(globalCellCells, adjncy, xadj);

    // decomposition options. 0 = use defaults
    List<int> options(3, 0);
    //options[0] = 1;     // don't use defaults but use values below
    //options[1] = -1;    // full debug info
    //options[2] = 15;    // random number seed

    // cell weights (so on the vertices of the dual)
    Field<int> cellWeights;

    // face weights (so on the edges of the dual)
    Field<int> faceWeights;

    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("metisCoeffs"))
    {
        const dictionary& metisCoeffs =
            decompositionDict_.subDict("metisCoeffs");
        word weightsFile;

        if (metisCoeffs.readIfPresent("cellWeightsFile", weightsFile))
        {
            Info<< "parMetisDecomp : Using cell-based weights read from "
                << weightsFile << endl;

            labelIOField cellIOWeights
            (
                IOobject
                (
                    weightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
            cellWeights.transfer(cellIOWeights);

            if (cellWeights.size() != cellCentres.size())
            {
                FatalErrorIn
                (
                    "parMetisDecomp::decompose"
                    "(const labelListList&, const pointField&)"
                )   << "Number of cell weights " << cellWeights.size()
                    << " read from " << cellIOWeights.objectPath()
                    << " does not equal number of cells " << cellCentres.size()
                    << exit(FatalError);
            }
        }

        //- faceWeights disabled. Only makes sense for cellCells from mesh.
        //if (metisCoeffs.readIfPresent("faceWeightsFile", weightsFile))
        //{
        //    Info<< "parMetisDecomp : Using face-based weights read from "
        //        << weightsFile << endl;
        //
        //    labelIOField weights
        //    (
        //        IOobject
        //        (
        //            weightsFile,
        //            mesh_.time().timeName(),
        //            mesh_,
        //            IOobject::MUST_READ,
        //            IOobject::AUTO_WRITE
        //        )
        //    );
        //
        //    if (weights.size() != mesh_.nFaces())
        //    {
        //        FatalErrorIn("parMetisDecomp::decompose(const pointField&)")
        //            << "Number of face weights " << weights.size()
        //            << " does not equal number of internal and boundary faces "
        //            << mesh_.nFaces()
        //            << exit(FatalError);
        //    }
        //
        //    faceWeights.setSize(2*mesh_.nInternalFaces()+nCoupledFaces);
        //
        //    // Assume symmetric weights. Keep same ordering as adjncy.
        //    nFacesPerCell = 0;
        //
        //    // Handle internal faces
        //    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        //    {
        //        label w = weights[faceI];
        //
        //        label own = faceOwner[faceI];
        //        label nei = faceNeighbour[faceI];
        //
        //        faceWeights[xadj[own] + nFacesPerCell[own]++] = w;
        //        faceWeights[xadj[nei] + nFacesPerCell[nei]++] = w;
        //    }
        //    // Coupled boundary faces
        //    forAll(patches, patchI)
        //    {
        //        const polyPatch& pp = patches[patchI];
        //
        //        if (pp.coupled())
        //        {
        //            label faceI = pp.start();
        //
        //            forAll(pp, i)
        //            {
        //                label w = weights[faceI];
        //                label own = faceOwner[faceI];
        //                adjncy[xadj[own] + nFacesPerCell[own]++] = w;
        //                faceI++;
        //            }
        //        }
        //    }
        //}

        if (metisCoeffs.readIfPresent("options", options))
        {
            Info<< "Using Metis options     " << options
                << nl << endl;

            if (options.size() != 3)
            {
                FatalErrorIn
                (
                    "parMetisDecomp::decompose"
                    "(const labelListList&, const pointField&)"
                )   << "Number of options " << options.size()
                    << " should be three." << exit(FatalError);
            }
        }
    }


    // Do actual decomposition
    List<int> finalDecomp;
    decompose
    (
        xadj,
        adjncy,
        cellCentres,
        cellWeights,
        faceWeights,
        options,

        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


// ************************************************************************* //
