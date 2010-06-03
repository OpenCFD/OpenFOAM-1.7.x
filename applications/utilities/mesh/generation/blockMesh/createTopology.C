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

\*---------------------------------------------------------------------------*/

#include "blockMesh.H"
#include "Time.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::blockMesh::blockLabelsOK
(
    const label blockLabel,
    const pointField& points,
    const cellShape& blockShape
)
{
    bool ok = true;

    forAll(blockShape, blockI)
    {
        if (blockShape[blockI] < 0)
        {
            ok = false;

            WarningIn
            (
                "bool Foam::blockMesh::blockLabelsOK"
                "(const label blockLabel, const pointField& points, "
                "const cellShape& blockShape)"
            )   << "block " << blockLabel
                << " point label " << blockShape[blockI]
                << " less than zero" << endl;
        }
        else if (blockShape[blockI] >= points.size())
        {
            ok = false;

            WarningIn
            (
                "bool Foam::blockMesh::blockLabelsOK"
                "(const label blockLabel, const pointField& points, "
                "const cellShape& blockShape)"
            )   << "block " << blockLabel
                << " point label " << blockShape[blockI]
                << " larger than " << points.size() - 1
                << " the largest defined point label" << endl;
        }
    }

    return ok;
}


bool Foam::blockMesh::patchLabelsOK
(
    const label patchLabel,
    const pointField& points,
    const faceList& patchFaces
)
{
    bool ok = true;

    forAll(patchFaces, faceI)
    {
        const labelList& f = patchFaces[faceI];

        forAll(f, fp)
        {
            if (f[fp] < 0)
            {
                ok = false;

                WarningIn
                (
                    "bool Foam::blockMesh::patchLabelsOK(...)"
                )   << "patch " << patchLabel
                    << " face " << faceI
                    << " point label " << f[fp]
                    << " less than zero" << endl;
            }
            else if (f[fp] >= points.size())
            {
                ok = false;

                WarningIn
                (
                    "bool Foam::blockMesh::patchLabelsOK(...)"
                )   << "patch " << patchLabel
                    << " face " << faceI
                    << " point label " << f[fp]
                    << " larger than " << points.size() - 1
                    << " the largest defined point label" << endl;
            }
        }
    }

    return ok;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::polyMesh* Foam::blockMesh::createTopology(IOdictionary& meshDescription)
{
    bool topologyOK = true;

    blockMesh& blocks = *this;

    word defaultPatchName = "defaultFaces";
    word defaultPatchType = emptyPolyPatch::typeName;

    // get names/types for the unassigned patch faces
    // this is a bit heavy handed (and ugly), but there is currently
    // no easy way to rename polyMesh patches subsequently
    if (const dictionary* dictPtr = meshDescription.subDictPtr("defaultPatch"))
    {
        dictPtr->readIfPresent("name", defaultPatchName);
        dictPtr->readIfPresent("type", defaultPatchType);
    }

    Info<< nl << "Creating blockCorners" << endl;

    // create blockCorners
    pointField tmpBlockPoints(meshDescription.lookup("vertices"));

    if (meshDescription.found("edges"))
    {
        // read number of non-linear edges in mesh
        Info<< nl << "Creating curved edges" << endl;

        ITstream& edgesStream(meshDescription.lookup("edges"));

        label nEdges = 0;

        token firstToken(edgesStream);

        if (firstToken.isLabel())
        {
            nEdges = firstToken.labelToken();
            edges_.setSize(nEdges);
        }
        else
        {
            edgesStream.putBack(firstToken);
        }

        // Read beginning of edges
        edgesStream.readBegin("edges");

        nEdges = 0;

        token lastToken(edgesStream);
        while
        (
            !(
                lastToken.isPunctuation()
                && lastToken.pToken() == token::END_LIST
            )
        )
        {
            if (edges_.size() <= nEdges)
            {
                edges_.setSize(nEdges + 1);
            }

            edgesStream.putBack(lastToken);

            edges_.set
            (
                nEdges,
                curvedEdge::New(tmpBlockPoints, edgesStream)
            );

            nEdges++;

            edgesStream >> lastToken;
        }
        edgesStream.putBack(lastToken);

        // Read end of edges
        edgesStream.readEnd("edges");
    }
    else
    {
        Info<< nl << "There are no non-linear edges" << endl;
    }


    Info<< nl << "Creating blocks" << endl;
    {
        ITstream& blockDescriptorStream(meshDescription.lookup("blocks"));

        // read number of blocks in mesh
        label nBlocks = 0;

        token firstToken(blockDescriptorStream);

        if (firstToken.isLabel())
        {
            nBlocks = firstToken.labelToken();
            blocks.setSize(nBlocks);
        }
        else
        {
            blockDescriptorStream.putBack(firstToken);
        }

        // Read beginning of blocks
        blockDescriptorStream.readBegin("blocks");

        nBlocks = 0;

        token lastToken(blockDescriptorStream);
        while
        (
            !(
                lastToken.isPunctuation()
                && lastToken.pToken() == token::END_LIST
            )
        )
        {
            if (blocks.size() <= nBlocks)
            {
                blocks.setSize(nBlocks + 1);
            }

            blockDescriptorStream.putBack(lastToken);

            blocks.set
            (
                nBlocks,
                new block
                (
                    blockDescriptor
                    (
                        tmpBlockPoints,
                        edges_,
                        blockDescriptorStream
                    )
                )
            );

            topologyOK = topologyOK && blockLabelsOK
            (
                nBlocks,
                tmpBlockPoints,
                blocks[nBlocks].blockDef().blockShape()
            );

            nBlocks++;

            blockDescriptorStream >> lastToken;
        }
        blockDescriptorStream.putBack(lastToken);

        // Read end of blocks
        blockDescriptorStream.readEnd("blocks");
    }


    Info<< nl << "Creating patches" << endl;

    faceListList tmpBlocksPatches;
    wordList patchNames;
    wordList patchTypes;

    {
        ITstream& patchStream(meshDescription.lookup("patches"));

        // read number of patches in mesh
        label nPatches = 0;

        token firstToken(patchStream);

        if (firstToken.isLabel())
        {
            nPatches = firstToken.labelToken();

            tmpBlocksPatches.setSize(nPatches);
            patchNames.setSize(nPatches);
            patchTypes.setSize(nPatches);
        }
        else
        {
            patchStream.putBack(firstToken);
        }

        // Read beginning of blocks
        patchStream.readBegin("patches");

        nPatches = 0;

        token lastToken(patchStream);
        while
        (
            !(
                lastToken.isPunctuation()
                && lastToken.pToken() == token::END_LIST
            )
        )
        {
            if (tmpBlocksPatches.size() <= nPatches)
            {
                tmpBlocksPatches.setSize(nPatches + 1);
                patchNames.setSize(nPatches + 1);
                patchTypes.setSize(nPatches + 1);
            }

            patchStream.putBack(lastToken);

            patchStream
                >> patchTypes[nPatches]
                >> patchNames[nPatches]
                >> tmpBlocksPatches[nPatches];


            // Catch multiple patches asap.
            for (label i = 0; i < nPatches; i++)
            {
                if (patchNames[nPatches] == patchNames[i])
                {
                    FatalErrorIn
                    (
                        "blockMesh::createTopology(IOdictionary&)"
                    )   << "Duplicate patch " << patchNames[nPatches]
                        << " at line " << patchStream.lineNumber()
                        << ". Exiting !" << nl
                        << exit(FatalError);
                }
            }

            topologyOK = topologyOK && patchLabelsOK
            (
                nPatches,
                tmpBlockPoints,
                tmpBlocksPatches[nPatches]
            );

            nPatches++;

            patchStream >> lastToken;
        }
        patchStream.putBack(lastToken);

        // Read end of blocks
        patchStream.readEnd("patches");
    }


    if (!topologyOK)
    {
        FatalErrorIn("blockMesh::createTopology(IOdictionary&)")
            << "Cannot create mesh due to errors in topology, exiting !" << nl
            << exit(FatalError);
    }


    Info<< nl << "Creating block mesh topology" << endl;

    PtrList<cellShape> tmpBlockCells(blocks.size());
    forAll(blocks, blockLabel)
    {
        tmpBlockCells.set
        (
            blockLabel,
            new cellShape(blocks[blockLabel].blockDef().blockShape())
        );

        if (tmpBlockCells[blockLabel].mag(tmpBlockPoints) < 0.0)
        {
            WarningIn
            (
                "blockMesh::createTopology(IOdictionary&)"
            )   << "negative volume block : " << blockLabel
                << ", probably defined inside-out" << endl;
        }
    }

    wordList patchPhysicalTypes(tmpBlocksPatches.size());

    preservePatchTypes
    (
        meshDescription.time(),
        meshDescription.time().constant(),
        polyMesh::meshSubDir,
        patchNames,
        patchTypes,
        defaultPatchName,
        defaultPatchType,
        patchPhysicalTypes
    );

    polyMesh* blockMeshPtr = new polyMesh
    (
        IOobject
        (
            "blockMesh",
            meshDescription.time().constant(),
            meshDescription.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        xferMove(tmpBlockPoints),
        tmpBlockCells,
        tmpBlocksPatches,
        patchNames,
        patchTypes,
        defaultPatchName,
        defaultPatchType,
        patchPhysicalTypes
    );

    checkBlockMesh(*blockMeshPtr);

    return blockMeshPtr;
}


// ************************************************************************* //
