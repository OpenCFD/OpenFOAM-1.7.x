/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

    From scotch forum:

    By: Francois PELLEGRINI RE: Graph mapping 'strategy' string [ reply ]
    2008-08-22 10:09 Strategy handling in Scotch is a bit tricky. In order
    not to be confused, you must have a clear view of how they are built.
    Here are some rules:

    1- Strategies are made up of "methods" which are combined by means of
    "operators".

    2- A method is of the form "m{param=value,param=value,...}", where "m"
    is a single character (this is your first error: "f" is a method name,
    not a parameter name).

    3- There exist different sort of strategies : bipartitioning strategies,
    mapping strategies, ordering strategies, which cannot be mixed. For
    instance, you cannot build a bipartitioning strategy and feed it to a
    mapping method (this is your second error).

    To use the "mapCompute" routine, you must create a mapping strategy, not
    a bipartitioning one, and so use stratGraphMap() and not
    stratGraphBipart(). Your mapping strategy should however be based on the
    "recursive bipartitioning" method ("b"). For instance, a simple (and
    hence not very efficient) mapping strategy can be :

    "b{sep=f}"

    which computes mappings with the recursive bipartitioning method "b",
    this latter using the Fiduccia-Mattheyses method "f" to compute its
    separators.

    If you want an exact partition (see your previous post), try
    "b{sep=fx}".

    However, these strategies are not the most efficient, as they do not
    make use of the multi-level framework.

    To use the multi-level framework, try for instance:

    "b{sep=m{vert=100,low=h,asc=f}x}"

    The current default mapping strategy in Scotch can be seen by using the
    "-vs" option of program gmap. It is, to date:

    b
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
                    org=f{move=80,pass=-1,bal=0.005},
                    width=3
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
                    width=3
                },
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
#include "Time.H"

static const char* notImplementedMessage =
"You are trying to use scotch but do not have the scotchDecomp library loaded."
"\nThis message is from the dummy scotchDecomp stub library instead.\n"
"\n"
"Please install scotch and make sure that libscotch.so is in your "
"LD_LIBRARY_PATH.\n"
"The scotchDecomp library can then be built in $FOAM_SRC/decompositionMethods/"
"scotchDecomp\n";


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
{}


Foam::label Foam::scotchDecomp::decompose
(
    const List<int>& adjncy,
    const List<int>& xadj,
    const scalarField& cWeights,

    List<int>& finalDecomp
)
{
    FatalErrorIn
    (
        "label scotchDecomp::decompose"
        "("
            "const List<int>&, "
            "const List<int>&, "
            "const scalarField&, "
            "List<int>&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return -1;
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
    FatalErrorIn
    (
        "labelList scotchDecomp::decompose"
        "("
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
)
{
    FatalErrorIn
    (
        "labelList scotchDecomp::decompose"
        "("
            "const labelList&, "
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    FatalErrorIn
    (
        "labelList scotchDecomp::decompose"
        "("
            "const labelListList&, "
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


void Foam::scotchDecomp::calcCSR
(
    const polyMesh& mesh,
    List<int>& adjncy,
    List<int>& xadj
)
{
    FatalErrorIn
    (
        "labelList scotchDecomp::decompose"
        "("
            "const polyMesh&, "
            "const List<int>&, "
            "const List<int>&"
        ")"
    )   << notImplementedMessage << exit(FatalError);
}


void Foam::scotchDecomp::calcCSR
(
    const labelListList& cellCells,
    List<int>& adjncy,
    List<int>& xadj
)
{
    FatalErrorIn
    (
        "labelList scotchDecomp::decompose"
        "("
            "const labelListList&, "
            "const List<int>&, "
            "const List<int>&"
        ")"
    )   << notImplementedMessage << exit(FatalError);
}


// ************************************************************************* //
