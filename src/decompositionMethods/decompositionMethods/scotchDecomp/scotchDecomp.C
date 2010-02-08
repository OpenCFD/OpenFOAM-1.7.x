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
#include "Time.H"
#include "cyclicPolyPatch.H"
#include "OFstream.H"
#include "metisDecomp.H"

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
            << "Call to scotch routine " << str << " failed."
            << exit(FatalError);
    }
}


// Call scotch with options from dictionary.
Foam::label Foam::scotchDecomp::decompose
(
    const List<int>& adjncy,
    const List<int>& xadj,
    const scalarField& cWeights,

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

    List<int> velotab;


    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "scotchDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != xadj.size()-1)
        {
            FatalErrorIn
            (
                "scotchDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << xadj.size()-1
                << exit(FatalError);
        }

        // Convert to integers.
        velotab.setSize(cWeights.size());
        forAll(velotab, i)
        {
            velotab[i] = int(cWeights[i]/minWeights);
        }
    }



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
            velotab.begin(),        // velotab, vertex weights
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

Foam::labelList Foam::scotchDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "scotchDecomp::decompose(const pointField&, const scalarField&)"
        )
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
    List<int> adjncy;
    List<int> xadj;
    metisDecomp::calcMetisCSR
    (
        mesh_,
        adjncy,
        xadj
    );

    // Decompose using default weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, pointWeights, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
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

        metisDecomp::calcMetisCSR(cellCells, adjncy, xadj);
    }

    // Decompose using weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, pointWeights, finalDecomp);

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
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "scotchDecomp::decompose"
            "(const labelListList&, const pointField&, const scalarField&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    List<int> adjncy;
    List<int> xadj;
    metisDecomp::calcMetisCSR(globalCellCells, adjncy, xadj);

    // Decompose using weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, cWeights, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


// ************************************************************************* //
