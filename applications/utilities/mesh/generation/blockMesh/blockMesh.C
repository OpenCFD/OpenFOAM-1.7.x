/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    blockMesh

Description
    Mesh generator

\*---------------------------------------------------------------------------*/

#include "blockMesh.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOdictionary
Foam::blockMesh::blockMesh(IOdictionary& meshDescription)
:
    topologyPtr_(createTopology(meshDescription)),
    blockOffsets_(createBlockOffsets()),
    mergeList_(createMergeList()),
    points_(createPoints(meshDescription)),
    cells_(createCells()),
    patches_(createPatches())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockMesh::~blockMesh()
{
    delete topologyPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::blockMesh::topology() const
{
    if (!topologyPtr_)
    {
        FatalErrorIn("blockMesh::topology() const")
            << "topologyPtr_ not allocated"
            << exit(FatalError);
    }

    return *topologyPtr_;
}


Foam::wordList Foam::blockMesh::patchNames() const
{
    const polyPatchList& patchTopologies = topology().boundaryMesh();
    wordList names(patchTopologies.size());

    forAll (names, patchI)
    {
        names[patchI] = patchTopologies[patchI].name();
    }

    return names;
}


Foam::wordList Foam::blockMesh::patchTypes() const
{
    const polyPatchList& patchTopologies = topology().boundaryMesh();
    wordList types(patchTopologies.size());

    forAll (types, patchI)
    {
        types[patchI] = patchTopologies[patchI].type();
    }

    return types;
}


Foam::wordList Foam::blockMesh::patchPhysicalTypes() const
{
    const polyPatchList& patchTopologies = topology().boundaryMesh();
    wordList physicalTypes(patchTopologies.size());

    forAll (physicalTypes, patchI)
    {
        physicalTypes[patchI] = patchTopologies[patchI].physicalType();
    }

    return physicalTypes;
}


Foam::label Foam::blockMesh::numZonedBlocks() const
{
    label num = 0;

    forAll(*this, blockI)
    {
        if (operator[](blockI).blockDef().zoneName().size())
        {
            num++;
        }
    }

    return num;
}


void Foam::blockMesh::writeTopology(Ostream& os) const
{
    const pointField& pts = topology().points();

    forAll(pts, pI)
    {
        const point& pt = pts[pI];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }

    const edgeList& edges = topology().edges();

    forAll(edges, eI)
    {
        const edge& e = edges[eI];

        os << "l " << e.start() + 1 << ' ' << e.end() + 1 << endl;
    }
}

// ************************************************************************* //
