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

#include "parMetisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "Time.H"

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

static const char* notImplementedMessage =
"You are trying to use parMetis but do not have the parMetisDecomp library "
"loaded.\n"
"This message is from the dummy parMetisDecomp stub library instead.\n"        
"\n"                                                                           
"Please install parMetis and make sure that libparMetis.so is in your "            
"LD_LIBRARY_PATH.\n"                                                           
"The parMetisDecomp library can then be built in $FOAM_SRC/decompositionMethods/"
"parMetisDecomp\n";                                                              


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
    FatalErrorIn
    (
        "label parMetisDecomp::decompose"
        "("
            "Field<int>&, "
            "Field<int>&, "
            "const pointField&, "
            "Field<int>&, "
            "Field<int>&, "
            "const List<int>&, "
            "List<int>&"
        ")"
    )<< notImplementedMessage << exit(FatalError);

    return -1;
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

Foam::labelList Foam::parMetisDecomp::decompose
(
    const pointField& cc,
    const scalarField& cWeights
)
{
    FatalErrorIn
    (
        "labelList parMetisDecomp::decompose"
        "("
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList();
}


Foam::labelList Foam::parMetisDecomp::decompose
(
    const labelList& cellToRegion,
    const pointField& regionPoints,
    const scalarField& regionWeights
)
{
    FatalErrorIn
    (
        "labelList parMetisDecomp::decompose"
        "("
            "const labelList&, "
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList();
}


Foam::labelList Foam::parMetisDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    FatalErrorIn
    (
        "labelList parMetisDecomp::decompose"
        "("
            "const labelListList&, "
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList();
}


void Foam::parMetisDecomp::calcMetisDistributedCSR
(
    const polyMesh& mesh,
    List<int>& adjncy,
    List<int>& xadj
)
{
    FatalErrorIn
    (
        "void parMetisDecomp::calcMetisDistributedCSR"
        "("
            "const polyMesh&, "
            "List<int>&, "
            "List<int>&"
        ")"
    )   << notImplementedMessage << exit(FatalError);
}


// ************************************************************************* //
