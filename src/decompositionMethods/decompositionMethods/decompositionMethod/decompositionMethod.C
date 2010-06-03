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

InClass
    decompositionMethod

\*---------------------------------------------------------------------------*/

#include "decompositionMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionMethod, 0);
    defineRunTimeSelectionTable(decompositionMethod, dictionary);
    defineRunTimeSelectionTable(decompositionMethod, dictionaryMesh);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
    const dictionary& decompositionDict
)
{
    word decompositionMethodTypeName(decompositionDict.lookup("method"));

    Info<< "Selecting decompositionMethod "
        << decompositionMethodTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(decompositionMethodTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "decompositionMethod::New"
            "(const dictionary& decompositionDict)"
        )   << "Unknown decompositionMethod "
            << decompositionMethodTypeName << endl << endl
            << "Valid decompositionMethods are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<decompositionMethod>(cstrIter()(decompositionDict));
}


Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
{
    word decompositionMethodTypeName(decompositionDict.lookup("method"));

    Info<< "Selecting decompositionMethod "
        << decompositionMethodTypeName << endl;

    dictionaryMeshConstructorTable::iterator cstrIter =
        dictionaryMeshConstructorTablePtr_->find(decompositionMethodTypeName);

    if (cstrIter == dictionaryMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "decompositionMethod::New"
            "(const dictionary& decompositionDict, "
            "const polyMesh& mesh)"
        )   << "Unknown decompositionMethod "
            << decompositionMethodTypeName << endl << endl
            << "Valid decompositionMethods are : " << endl
            << dictionaryMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<decompositionMethod>(cstrIter()(decompositionDict, mesh));
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const pointField& points
)
{
    scalarField weights(0);

    return decompose(points, weights);
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const labelList& fineToCoarse,
    const pointField& coarsePoints,
    const scalarField& coarseWeights
)
{
    // Decompose based on agglomerated points
    labelList coarseDistribution(decompose(coarsePoints, coarseWeights));

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const labelList& fineToCoarse,
    const pointField& coarsePoints
)
{
    // Decompose based on agglomerated points
    labelList coarseDistribution(decompose(coarsePoints));

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
    }

    return fineDistribution;
}


void Foam::decompositionMethod::calcCellCells
(
    const polyMesh& mesh,
    const labelList& fineToCoarse,
    const label nCoarse,
    labelListList& cellCells
)
{
    if (fineToCoarse.size() != mesh.nCells())
    {
        FatalErrorIn
        (
            "decompositionMethod::calcCellCells"
            "(const labelList&, labelListList&) const"
        )   << "Only valid for mesh agglomeration." << exit(FatalError);
    }

    List<DynamicList<label> > dynCellCells(nCoarse);

    forAll(mesh.faceNeighbour(), faceI)
    {
        label own = fineToCoarse[mesh.faceOwner()[faceI]];
        label nei = fineToCoarse[mesh.faceNeighbour()[faceI]];

        if (own != nei)
        {
            if (findIndex(dynCellCells[own], nei) == -1)
            {
                dynCellCells[own].append(nei);
            }
            if (findIndex(dynCellCells[nei], own) == -1)
            {
                dynCellCells[nei].append(own);
            }
        }
    }

    cellCells.setSize(dynCellCells.size());
    forAll(dynCellCells, coarseI)
    {
        cellCells[coarseI].transfer(dynCellCells[coarseI]);
    }
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const labelListList& globalCellCells,
    const pointField& cc
)
{
    scalarField cWeights(0);

    return decompose(globalCellCells, cc, cWeights);
}


// ************************************************************************* //
