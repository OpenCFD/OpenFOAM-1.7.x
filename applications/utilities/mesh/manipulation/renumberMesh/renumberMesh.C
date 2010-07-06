/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anispulation  |
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

Application
    renumberMesh

Description
    Renumbers the cell list in order to reduce the bandwidth, reading and
    renumbering all fields from all the time directories.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "polyTopoChange.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "bandCompression.H"
#include "faceSet.H"
#include "SortableList.H"
#include "decompositionMethod.H"
#include "fvMeshSubset.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam;


// Calculate band of matrix
label getBand(const labelList& owner, const labelList& neighbour)
{
    label band = 0;

    forAll(neighbour, faceI)
    {
        label diff = neighbour[faceI] - owner[faceI];

        if (diff > band)
        {
            band = diff;
        }
    }
    return band;
}


// Return new to old cell numbering
labelList regionBandCompression
(
    const fvMesh& mesh,
    const labelList& cellToRegion
)
{
    Pout<< "Determining cell order:" << endl;

    labelList cellOrder(cellToRegion.size());

    label nRegions = max(cellToRegion)+1;

    labelListList regionToCells(invertOneToMany(nRegions, cellToRegion));

    label cellI = 0;

    forAll(regionToCells, regionI)
    {
        Pout<< "    region " << regionI << " starts at " << cellI << endl;

        // Per region do a reordering.
        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(cellToRegion, regionI);
        const fvMesh& subMesh = subsetter.subMesh();
        labelList subCellOrder(bandCompression(subMesh.cellCells()));

        const labelList& cellMap = subsetter.cellMap();

        forAll(subCellOrder, i)
        {
            cellOrder[cellI++] = cellMap[subCellOrder[i]];
        }
    }
    Pout<< endl;

    return cellOrder;
}


// Determine face order such that inside region faces are sorted
// upper-triangular but inbetween region faces are handled like boundary faces.
labelList regionFaceOrder
(
    const primitiveMesh& mesh,
    const labelList& cellOrder,     // new to old cell
    const labelList& cellToRegion   // old cell to region
)
{
    Pout<< "Determining face order:" << endl;

    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));

    labelList oldToNewFace(mesh.nFaces(), -1);

    label newFaceI = 0;

    label prevRegion = -1;

    forAll (cellOrder, newCellI)
    {
        label oldCellI = cellOrder[newCellI];

        if (cellToRegion[oldCellI] != prevRegion)
        {
            prevRegion = cellToRegion[oldCellI];
            Pout<< "    region " << prevRegion << " internal faces start at "
                << newFaceI << endl;
        }

        const cell& cFaces = mesh.cells()[oldCellI];

        SortableList<label> nbr(cFaces.size());

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if (mesh.isInternalFace(faceI))
            {
                // Internal face. Get cell on other side.
                label nbrCellI = reverseCellOrder[mesh.faceNeighbour()[faceI]];
                if (nbrCellI == newCellI)
                {
                    nbrCellI = reverseCellOrder[mesh.faceOwner()[faceI]];
                }

                if (cellToRegion[oldCellI] != cellToRegion[cellOrder[nbrCellI]])
                {
                    // Treat like external face. Do later.
                    nbr[i] = -1;
                }
                else if (newCellI < nbrCellI)
                {
                    // CellI is master
                    nbr[i] = nbrCellI;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNewFace[cFaces[nbr.indices()[i]]] = newFaceI++;
            }
        }
    }

    // Do region interfaces
    label nRegions = max(cellToRegion)+1;
    {
        // Sort in increasing region
        SortableList<label> sortKey(mesh.nFaces(), labelMax);

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            label ownRegion = cellToRegion[mesh.faceOwner()[faceI]];
            label neiRegion = cellToRegion[mesh.faceNeighbour()[faceI]];

            if (ownRegion != neiRegion)
            {
                sortKey[faceI] =
                    min(ownRegion, neiRegion)*nRegions
                   +max(ownRegion, neiRegion);
            }
        }
        sortKey.sort();

        // Extract.
        label prevKey = -1;
        forAll(sortKey, i)
        {
            label key = sortKey[i];

            if (key == labelMax)
            {
                break;
            }

            if (prevKey != key)
            {
                Pout<< "    faces inbetween region " << key/nRegions
                    << " and " << key%nRegions
                    << " start at " << newFaceI << endl;
                prevKey = key;
            }

            oldToNewFace[sortKey.indices()[i]] = newFaceI++;
        }
    }

    // Leave patch faces intact.
    for (label faceI = newFaceI; faceI < mesh.nFaces(); faceI++)
    {
        oldToNewFace[faceI] = faceI;
    }


    // Check done all faces.
    forAll(oldToNewFace, faceI)
    {
        if (oldToNewFace[faceI] == -1)
        {
            FatalErrorIn
            (
                "polyDualMesh::getFaceOrder"
                "(const labelList&, const labelList&, const label) const"
            )   << "Did not determine new position"
                << " for face " << faceI
                << abort(FatalError);
        }
    }
    Pout<< endl;

    return invert(mesh.nFaces(), oldToNewFace);
}


// cellOrder: old cell for every new cell
// faceOrder: old face for every new face. Ordering of boundary faces not
// changed.
autoPtr<mapPolyMesh> reorderMesh
(
    polyMesh& mesh,
    const labelList& cellOrder,
    const labelList& faceOrder
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));
    labelList reverseFaceOrder(invert(faceOrder.size(), faceOrder));

    faceList newFaces(reorder(reverseFaceOrder, mesh.faces()));
    labelList newOwner
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceOwner())
        )
    );
    labelList newNeighbour
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceNeighbour())
        )
    );

    // Check if any faces need swapping.
    forAll(newNeighbour, faceI)
    {
        label own = newOwner[faceI];
        label nei = newNeighbour[faceI];

        if (nei < own)
        {
            newFaces[faceI] = newFaces[faceI].reverseFace();
            Swap(newOwner[faceI], newNeighbour[faceI]);
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelList patchSizes(patches.size());
    labelList patchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    labelListList patchPointMap(patches.size());
    forAll(patches, patchI)
    {
        patchSizes[patchI] = patches[patchI].size();
        patchStarts[patchI] = patches[patchI].start();
        oldPatchNMeshPoints[patchI] = patches[patchI].nPoints();
        patchPointMap[patchI] = identity(patches[patchI].nPoints());
    }

    mesh.resetPrimitives
    (
        Xfer<pointField>::null(),
        xferMove(newFaces),
        xferMove(newOwner),
        xferMove(newNeighbour),
        patchSizes,
        patchStarts
    );

    return autoPtr<mapPolyMesh>
    (
        new mapPolyMesh
        (
            mesh,                       //const polyMesh& mesh,
            mesh.nPoints(),             // nOldPoints,
            mesh.nFaces(),              // nOldFaces,
            mesh.nCells(),              // nOldCells,
            identity(mesh.nPoints()),   // pointMap,
            List<objectMap>(0),         // pointsFromPoints,
            faceOrder,                  // faceMap,
            List<objectMap>(0),         // facesFromPoints,
            List<objectMap>(0),         // facesFromEdges,
            List<objectMap>(0),         // facesFromFaces,
            cellOrder,                  // cellMap,
            List<objectMap>(0),         // cellsFromPoints,
            List<objectMap>(0),         // cellsFromEdges,
            List<objectMap>(0),         // cellsFromFaces,
            List<objectMap>(0),         // cellsFromCells,
            identity(mesh.nPoints()),   // reversePointMap,
            reverseFaceOrder,           // reverseFaceMap,
            reverseCellOrder,           // reverseCellMap,
            labelHashSet(0),            // flipFaceFlux,
            patchPointMap,              // patchPointMap,
            labelListList(0),           // pointZoneMap,
            labelListList(0),           // faceZonePointMap,
            labelListList(0),           // faceZoneFaceMap,
            labelListList(0),           // cellZoneMap,
            pointField(0),              // preMotionPoints,
            patchStarts,                // oldPatchStarts,
            oldPatchNMeshPoints         // oldPatchNMeshPoints
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("blockOrder", "");
    argList::validOptions.insert("orderPoints", "");
    argList::validOptions.insert("writeMaps", "");
    argList::validOptions.insert("overwrite", "");

#   include "addTimeOptions.H"

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool blockOrder = args.optionFound("blockOrder");
    if (blockOrder)
    {
        Info<< "Ordering cells into regions (using decomposition);"
            << " ordering faces into region-internal and region-external." << nl
            << endl;
    }

    const bool orderPoints = args.optionFound("orderPoints");
    if (orderPoints)
    {
        Info<< "Ordering points into internal and boundary points." << nl
            << endl;
    }

    const bool writeMaps = args.optionFound("writeMaps");

    if (writeMaps)
    {
        Info<< "Writing renumber maps (new to old) to polyMesh." << nl
            << endl;
    }

    bool overwrite = args.optionFound("overwrite");

    label band = getBand(mesh.faceOwner(), mesh.faceNeighbour());

    Info<< "Mesh size: " << returnReduce(mesh.nCells(), sumOp<label>()) << nl
        << "Band before renumbering: "
        << returnReduce(band, maxOp<label>()) << nl << endl;


    // Read parallel reconstruct maps
    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellProcAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        labelList(0)
    );

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceProcAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        labelList(0)
    );
    labelIOList pointProcAddressing
    (
        IOobject
        (
            "pointProcAddressing",
            mesh.pointsInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        labelList(0)
    );
    labelIOList boundaryProcAddressing
    (
        IOobject
        (
            "boundaryProcAddressing",
            mesh.pointsInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        labelList(0)
    );


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);


    autoPtr<mapPolyMesh> map;

    if (blockOrder)
    {
        // Renumbering in two phases. Should be done in one so mapping of
        // fields is done correctly!

        // Read decomposePar dictionary
        IOdictionary decomposeDict
        (
            IOobject
            (
                "decomposeParDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        autoPtr<decompositionMethod> decomposePtr = decompositionMethod::New
        (
            decomposeDict,
            mesh
        );

        labelList cellToRegion(decomposePtr().decompose(mesh.cellCentres()));

        // For debugging: write out region
        {
            volScalarField cellDist
            (
                IOobject
                (
                    "cellDist",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("cellDist", dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            );

            forAll(cellToRegion, cellI)
            {
               cellDist[cellI] = cellToRegion[cellI];
            }

            cellDist.write();

            Info<< nl << "Written decomposition as volScalarField to "
                << cellDist.name() << " for use in postprocessing."
                << nl << endl;
        }

        // Use block based renumbering.
        //labelList cellOrder(bandCompression(mesh.cellCells()));
        labelList cellOrder(regionBandCompression(mesh, cellToRegion));

        // Determine new to old face order with new cell numbering
        labelList faceOrder
        (
            regionFaceOrder
            (
                mesh,
                cellOrder,
                cellToRegion
            )
        );

        if (!overwrite)
        {
            runTime++;
        }

        // Change the mesh.
        map = reorderMesh(mesh, cellOrder, faceOrder);
    }
    else
    {
        // Use built-in renumbering.

        polyTopoChange meshMod(mesh);

        if (!overwrite)
        {
            runTime++;
        }

        // Change the mesh.
        map = meshMod.changeMesh
        (
            mesh,
            false,      // inflate
            true,       // parallel sync
            true,       // cell ordering
            orderPoints // point ordering
        );
    }

    // Update fields
    mesh.updateMesh(map);

    // Update proc maps
    if (cellProcAddressing.headerOk())
    {
        Info<< "Renumbering processor cell decomposition map "
            << cellProcAddressing.name() << endl;

        cellProcAddressing = labelList
        (
            UIndirectList<label>(cellProcAddressing, map().cellMap())
        );
    }
    if (faceProcAddressing.headerOk())
    {
        Info<< "Renumbering processor face decomposition map "
            << faceProcAddressing.name() << endl;

        faceProcAddressing = labelList
        (
            UIndirectList<label>(faceProcAddressing, map().faceMap())
        );
    }
    if (pointProcAddressing.headerOk())
    {
        Info<< "Renumbering processor point decomposition map "
            << pointProcAddressing.name() << endl;

        pointProcAddressing = labelList
        (
            UIndirectList<label>(pointProcAddressing, map().pointMap())
        );
    }


    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }


    band = getBand(mesh.faceOwner(), mesh.faceNeighbour());

    Info<< "Band after renumbering: "
        << returnReduce(band, maxOp<label>()) << nl << endl;


    if (orderPoints)
    {
        // Force edge calculation (since only reason that points would need to
        // be sorted)
        (void)mesh.edges();

        label nTotPoints = returnReduce
        (
            mesh.nPoints(),
            sumOp<label>()
        );
        label nTotIntPoints = returnReduce
        (
            mesh.nInternalPoints(),
            sumOp<label>()
        );

        label nTotEdges = returnReduce
        (
            mesh.nEdges(),
            sumOp<label>()
        );
        label nTotIntEdges = returnReduce
        (
            mesh.nInternalEdges(),
            sumOp<label>()
        );
        label nTotInt0Edges = returnReduce
        (
            mesh.nInternal0Edges(),
            sumOp<label>()
        );
        label nTotInt1Edges = returnReduce
        (
            mesh.nInternal1Edges(),
            sumOp<label>()
        );

        Info<< "Points:" << nl
            << "    total   : " << nTotPoints << nl
            << "    internal: " << nTotIntPoints << nl
            << "    boundary: " << nTotPoints-nTotIntPoints << nl
            << "Edges:" << nl
            << "    total   : " << nTotEdges << nl
            << "    internal: " << nTotIntEdges << nl
            << "        internal using 0 boundary points: "
            << nTotInt0Edges << nl
            << "        internal using 1 boundary points: "
            << nTotInt1Edges-nTotInt0Edges << nl
            << "        internal using 2 boundary points: "
            << nTotIntEdges-nTotInt1Edges << nl
            << "    boundary: " << nTotEdges-nTotIntEdges << nl
            << endl;
    }


    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh to " << mesh.facesInstance() << endl;

    mesh.write();
    if (cellProcAddressing.headerOk())
    {
        cellProcAddressing.instance() = mesh.facesInstance();
        cellProcAddressing.write();
    }
    if (faceProcAddressing.headerOk())
    {
        faceProcAddressing.instance() = mesh.facesInstance();
        faceProcAddressing.write();
    }
    if (pointProcAddressing.headerOk())
    {
        pointProcAddressing.instance() = mesh.facesInstance();
        pointProcAddressing.write();
    }
    if (boundaryProcAddressing.headerOk())
    {
        boundaryProcAddressing.instance() = mesh.facesInstance();
        boundaryProcAddressing.write();
    }


    if (writeMaps)
    {
        labelIOList
        (
            IOobject
            (
                "cellMap",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            map().cellMap()
        ).write();
        labelIOList
        (
            IOobject
            (
                "faceMap",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            map().faceMap()
        ).write();
        labelIOList
        (
            IOobject
            (
                "pointMap",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            map().pointMap()
        ).write();
    }

    Info<< "\nEnd.\n" << endl;

    return 0;
}


// ************************************************************************* //
