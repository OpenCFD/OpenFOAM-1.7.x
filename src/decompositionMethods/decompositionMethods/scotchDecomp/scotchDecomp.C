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

    Default strategy:
        Strat=b
        {
            job=t,
            map=t,
            poli=S,
            sep=
            (
                m
                {
                    asc=b
                    {
                        bnd=d{pass=40,dif=1,rem=1}f{move=80,pass=-1,bal=0.005},
                        org=f{move=80,pass=-1,bal=0.005},width=3
                    },
                    low=h{pass=10}f{move=80,pass=-1,bal=0.0005},
                    type=h,
                    vert=80,
                    rat=0.8
                }
              | m
                {
                    asc=b
                    {
                        bnd=d{pass=40,dif=1,rem=1}f{move=80,pass=-1,bal=0.005},
                        org=f{move=80,pass=-1,bal=0.005},
                        width=3},
                        low=h{pass=10}f{move=80,pass=-1,bal=0.0005},
                        type=h,
                        vert=80,
                        rat=0.8
                }
            )
        }
\*---------------------------------------------------------------------------*/

#include "scotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "IFstream.H"
#include "Time.H"
#include "cyclicPolyPatch.H"
#include "OFstream.H"

extern "C"
{
#define OMPI_SKIP_MPICXX
#include "module.h"
#include "common.h"
#include "scotch.h"
}


// Hack: scotch generates floating point errors so need to switch of error
//       trapping!
#if defined(linux) || defined(linuxAMD64) || defined(linuxIA64)
#    define LINUX
#endif

#if defined(LINUX) && defined(__GNUC__)
#    define LINUX_GNUC
#endif

#ifdef LINUX_GNUC
#   ifndef __USE_GNU
#       define __USE_GNU
#   endif
#   include <fenv.h>
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scotchDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotchDecomp,
        dictionaryMesh
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scotchDecomp::check(const int retVal, const char* str)
{
    if (retVal)
    {
        FatalErrorIn("scotchDecomp::decompose(..)")
            << "Called to scotch routine " << str << " failed."
            << exit(FatalError);
    }
}


// Call Metis with options from dictionary.
Foam::label Foam::scotchDecomp::decompose
(
    const List<int>& adjncy,
    const List<int>& xadj,

    List<int>& finalDecomp
)
{
    // Dump graph
    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");

        if (scotchCoeffs.found("writeGraph"))
        {
            Switch writeGraph(scotchCoeffs.lookup("writeGraph"));

            if (writeGraph)
            {
                OFstream str(mesh_.time().path() / mesh_.name() + ".grf");

                Info<< "Dumping Scotch graph file to " << str.name() << endl
                    << "Use this in combination with gpart." << endl;

                label version = 0;
                str << version << nl;
                // Numer of vertices
                str << xadj.size()-1 << ' ' << adjncy.size() << nl;
                // Numbering starts from 0
                label baseval = 0;
                // Has weights?
                label hasEdgeWeights = 0;
                label hasVertexWeights = 0;
                label numericflag = 10*hasEdgeWeights+hasVertexWeights;
                str << baseval << ' ' << numericflag << nl;
                for (label cellI = 0; cellI < xadj.size()-1; cellI++)
                {
                    label start = xadj[cellI];
                    label end = xadj[cellI+1];
                    str << end-start;

                    for (label i = start; i < end; i++)
                    {
                        str << ' ' << adjncy[i];
                    }
                    str << nl;
                }
            }
        }
    }


    // Strategy
    // ~~~~~~~~

    // Default.
    SCOTCH_Strat stradat;
    check(SCOTCH_stratInit(&stradat), "SCOTCH_stratInit");

    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");


        string strategy;
        if (scotchCoeffs.readIfPresent("strategy", strategy))
        {
            if (debug)
            {
                Info<< "scotchDecomp : Using strategy " << strategy << endl;
            }
            SCOTCH_stratGraphMap(&stradat, strategy.c_str());
            //fprintf(stdout, "S\tStrat=");
            //SCOTCH_stratSave(&stradat, stdout);
            //fprintf(stdout, "\n");
        }
    }


    // Graph
    // ~~~~~

    SCOTCH_Graph grafdat;
    check(SCOTCH_graphInit(&grafdat), "SCOTCH_graphInit");
    check
    (
        SCOTCH_graphBuild
        (
            &grafdat,
            0,                      // baseval, c-style numbering
            xadj.size()-1,          // vertnbr, nCells
            xadj.begin(),           // verttab, start index per cell into adjncy
            &xadj[1],               // vendtab, end index  ,,
            NULL,                   // velotab, vertex weights
            NULL,                   // vlbltab
            adjncy.size(),          // edgenbr, number of arcs
            adjncy.begin(),         // edgetab
            NULL                    // edlotab, edge weights
        ),
        "SCOTCH_graphBuild"
    );
    check(SCOTCH_graphCheck(&grafdat), "SCOTCH_graphCheck");


    // Architecture
    // ~~~~~~~~~~~~
    // (fully connected network topology since using switch)

    SCOTCH_Arch archdat;
    check(SCOTCH_archInit(&archdat), "SCOTCH_archInit");

    List<label> processorWeights;
    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");

        scotchCoeffs.readIfPresent("processorWeights", processorWeights);
    }
    if (processorWeights.size())
    {
        if (debug)
        {
            Info<< "scotchDecomp : Using procesor weights " << processorWeights
                << endl;
        }
        check
        (
            SCOTCH_archCmpltw(&archdat, nProcessors_, processorWeights.begin()),
            "SCOTCH_archCmpltw"
        );
    }
    else
    {
        check
        (
            SCOTCH_archCmplt(&archdat, nProcessors_),
            "SCOTCH_archCmplt"
        );
    }


    //SCOTCH_Mapping mapdat;
    //SCOTCH_graphMapInit(&grafdat, &mapdat, &archdat, NULL);
    //SCOTCH_graphMapCompute(&grafdat, &mapdat, &stradat); /* Perform mapping */
    //SCOTCH_graphMapExit(&grafdat, &mapdat);


    // Hack:switch off fpu error trapping
#   ifdef LINUX_GNUC
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
#   endif

    finalDecomp.setSize(xadj.size()-1);
    finalDecomp = 0;
    check
    (
        SCOTCH_graphMap
        (
            &grafdat,
            &archdat,
            &stradat,           // const SCOTCH_Strat *
            finalDecomp.begin() // parttab
        ),
        "SCOTCH_graphMap"
    );

#   ifdef LINUX_GNUC
    feenableexcept(oldExcepts);
#   endif



    //finalDecomp.setSize(xadj.size()-1);
    //check
    //(
    //    SCOTCH_graphPart
    //    (
    //        &grafdat,
    //        nProcessors_,       // partnbr
    //        &stradat,           // const SCOTCH_Strat *
    //        finalDecomp.begin() // parttab
    //    ),
    //    "SCOTCH_graphPart"
    //);

    // Release storage for graph
    SCOTCH_graphExit(&grafdat);
    // Release storage for strategy
    SCOTCH_stratExit(&stradat);
    // Release storage for network topology
    SCOTCH_archExit(&archdat);

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scotchDecomp::scotchDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::scotchDecomp::decompose(const pointField& points)
{
    if (points.size() != mesh_.nCells())
    {
        FatalErrorIn("scotchDecomp::decompose(const pointField&)")
            << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    List<int> xadj(mesh_.nCells()+1);

    // Initialise the number of internal faces of the cells to twice the
    // number of internal faces
    label nInternalFaces = 2*mesh_.nInternalFaces();

    // Check the boundary for coupled patches and add to the number of
    // internal faces
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchi)
    {
        if (isA<cyclicPolyPatch>(pbm[patchi]))
        {
            nInternalFaces += pbm[patchi].size();
        }
    }

    // Create the adjncy array the size of the total number of internal and
    // coupled faces
    List<int> adjncy(nInternalFaces);

    // Fill in xadj
    // ~~~~~~~~~~~~
    label freeAdj = 0;

    for (label cellI = 0; cellI < mesh_.nCells(); cellI++)
    {
        xadj[cellI] = freeAdj;

        const labelList& cFaces = mesh_.cells()[cellI];

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if
            (
                mesh_.isInternalFace(faceI)
             || isA<cyclicPolyPatch>(pbm[pbm.whichPatch(faceI)])
            )
            {
                freeAdj++;
            }
        }
    }
    xadj[mesh_.nCells()] = freeAdj;


    // Fill in adjncy
    // ~~~~~~~~~~~~~~

    labelList nFacesPerCell(mesh_.nCells(), 0);

    // Internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label nei = mesh_.faceNeighbour()[faceI];

        adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
        adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
    }

    // Coupled faces. Only cyclics done.
    forAll(pbm, patchi)
    {
        if (isA<cyclicPolyPatch>(pbm[patchi]))
        {
            const unallocLabelList& faceCells = pbm[patchi].faceCells();

            label sizeby2 = faceCells.size()/2;

            for (label facei=0; facei<sizeby2; facei++)
            {
                label own = faceCells[facei];
                label nei = faceCells[facei + sizeby2];

                adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
                adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
            }
        }
    }

    // Decompose using default weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


// From cell-cell connections to Metis format (like CompactListList)
void Foam::scotchDecomp::calcMetisCSR
(
    const labelListList& cellCells,
    List<int>& adjncy,
    List<int>& xadj
)
{
    // Count number of internal faces
    label nConnections = 0;

    forAll(cellCells, coarseI)
    {
        nConnections += cellCells[coarseI].size();
    }

    // Create the adjncy array as twice the size of the total number of
    // internal faces
    adjncy.setSize(nConnections);

    xadj.setSize(cellCells.size()+1);


    // Fill in xadj
    // ~~~~~~~~~~~~
    label freeAdj = 0;

    forAll(cellCells, coarseI)
    {
        xadj[coarseI] = freeAdj;

        const labelList& cCells = cellCells[coarseI];

        forAll(cCells, i)
        {
            adjncy[freeAdj++] = cCells[i];
        }
    }
    xadj[cellCells.size()] = freeAdj;
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelList& agglom,
    const pointField& agglomPoints
)
{
    if (agglom.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "parMetisDecomp::decompose(const labelList&, const pointField&)"
        )   << "Size of cell-to-coarse map " << agglom.size()
            << " differs from number of cells in mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    List<int> adjncy;
    List<int> xadj;
    {
        // Get cellCells on coarse mesh.
        labelListList cellCells;
        calcCellCells
        (
            mesh_,
            agglom,
            agglomPoints.size(),
            cellCells
        );

        calcMetisCSR(cellCells, adjncy, xadj);
    }

    // Decompose using default weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, finalDecomp);


    // Rework back into decomposition for original mesh_
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = finalDecomp[agglom[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "scotchDecomp::decompose(const pointField&, const labelListList&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    List<int> adjncy;
    List<int> xadj;
    calcMetisCSR(globalCellCells, adjncy, xadj);


    // Decompose using default weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


// ************************************************************************* //
