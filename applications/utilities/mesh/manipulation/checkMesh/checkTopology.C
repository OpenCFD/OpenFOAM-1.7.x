#include "checkTopology.H"
#include "polyMesh.H"
#include "Time.H"
#include "regionSplit.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOmanip.H"

bool Foam::checkSync(const wordList& names)
{
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = names;
    Pstream::gatherList(allNames);

    bool hasError = false;

    for (label procI = 1; procI < allNames.size(); procI++)
    {
        if (allNames[procI] != allNames[0])
        {
            hasError = true;

            Info<< " ***Inconsistent zones across processors, "
                   "processor 0 has zones:" << allNames[0]
                << ", processor " << procI << " has zones:"
                << allNames[procI]
                << endl;
        }
    }
    return hasError;
}


Foam::label Foam::checkTopology
(
    const polyMesh& mesh,
    const bool allTopology,
    const bool allGeometry
)
{
    label noFailedChecks = 0;

    Info<< "Checking topology..." << endl;

    // Check if the boundary definition is unique
    mesh.boundaryMesh().checkDefinition(true);

    // Check if the boundary processor patches are correct
    mesh.boundaryMesh().checkParallelSync(true);

    // Check names of zones are equal
    if (checkSync(mesh.cellZones().names()))
    {
        noFailedChecks++;
    }
    if (checkSync(mesh.faceZones().names()))
    {
        noFailedChecks++;
    }
    if (checkSync(mesh.pointZones().names()))
    {
        noFailedChecks++;
    }

    // Check contents of faceZones consistent
    {
        forAll(mesh.faceZones(), zoneI)
        {
            if (mesh.faceZones()[zoneI].checkParallelSync(false))
            {
                Info<< " ***FaceZone " << mesh.faceZones()[zoneI].name()
                    << " is not correctly synchronised"
                    << " across coupled boundaries."
                    << " (coupled faces both"
                    << " present in set but with opposite flipmap)" << endl;
                noFailedChecks++;
            }
        }
    }

    {
        pointSet points(mesh, "unusedPoints", mesh.nPoints()/100);
        if (mesh.checkPoints(true, &points))
        {
            noFailedChecks++;

            label nPoints = returnReduce(points.size(), sumOp<label>());

            Info<< "  <<Writing " << nPoints
                << " unused points to set " << points.name() << endl;
            points.write();
        }
    }

    {
        faceSet faces(mesh, "upperTriangularFace", mesh.nFaces()/100);
        if (mesh.checkUpperTriangular(true, &faces))
        {
            noFailedChecks++;
        }

        label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            Info<< "  <<Writing " << nFaces
                << " unordered faces to set " << faces.name() << endl;
            faces.write();
        }
    }

    if (allTopology)
    {
        cellSet cells(mesh, "zipUpCells", mesh.nCells()/100);
        if (mesh.checkCellsZipUp(true, &cells))
        {
            noFailedChecks++;

            label nCells = returnReduce(cells.size(), sumOp<label>());

            Info<< "  <<Writing " << nCells
                << " cells with over used edges to set " << cells.name()
                << endl;
            cells.write();
        }
    }

    {
        faceSet faces(mesh, "outOfRangeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceVertices(true, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            Info<< "  <<Writing " << nFaces
                << " faces with out-of-range or duplicate vertices to set "
                << faces.name() << endl;
            faces.write();
        }
    }

    if (allTopology)
    {
        faceSet faces(mesh, "edgeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceFaces(true, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            Info<< "  <<Writing " << nFaces
                << " faces with incorrect edges to set " << faces.name()
                << endl;
            faces.write();
        }
    }

    if (allTopology)
    {
        labelList nInternalFaces(mesh.nCells(), 0);

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            nInternalFaces[mesh.faceOwner()[faceI]]++;
            nInternalFaces[mesh.faceNeighbour()[faceI]]++;
        }
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        forAll(patches, patchI)
        {
            if (patches[patchI].coupled())
            {
                const unallocLabelList& owners = patches[patchI].faceCells();

                forAll(owners, i)
                {
                    nInternalFaces[owners[i]]++;
                }
            }
        }

        faceSet oneCells(mesh, "oneInternalFaceCells", mesh.nCells()/100);
        faceSet twoCells(mesh, "twoInternalFacesCells", mesh.nCells()/100);

        forAll(nInternalFaces, cellI)
        {
            if (nInternalFaces[cellI] <= 1)
            {
                oneCells.insert(cellI);
            }
            else if (nInternalFaces[cellI] == 2)
            {
                twoCells.insert(cellI);
            }
        }

        label nOneCells = returnReduce(oneCells.size(), sumOp<label>());

        if (nOneCells > 0)
        {
            Info<< "  <<Writing " << nOneCells
                << " cells with with single non-boundary face to set "
                << oneCells.name()
                << endl;
            oneCells.write();
        }

        label nTwoCells = returnReduce(twoCells.size(), sumOp<label>());

        if (nTwoCells > 0)
        {
            Info<< "  <<Writing " << nTwoCells
                << " cells with with single non-boundary face to set "
                << twoCells.name()
                << endl;
            twoCells.write();
        }
    }

    {
        regionSplit rs(mesh);

        if (rs.nRegions() == 1)
        {
            Info<< "    Number of regions: " << rs.nRegions() << " (OK)."
                << endl;

        }
        else
        {
            Info<< "   *Number of regions: " << rs.nRegions() << endl;

            Info<< "    The mesh has multiple regions which are not connected "
                   "by any face." << endl
                << "  <<Writing region information to "
                << mesh.time().timeName()/"cellToRegion"
                << endl;

            labelIOList ctr
            (
                IOobject
                (
                    "cellToRegion",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rs
            );
            ctr.write();
        }
    }

    if (!Pstream::parRun())
    {
        Pout<< "\nChecking patch topology for multiply connected surfaces ..."
            << endl;

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Non-manifold points
        pointSet points
        (
            mesh,
            "nonManifoldPoints",
            mesh.nPoints()/100
        );

        Pout.setf(ios_base::left);

        Pout<< "    "
            << setw(20) << "Patch"
            << setw(9) << "Faces"
            << setw(9) << "Points"
            << setw(34) << "Surface topology";
        if (allGeometry)
        {
            Pout<< " Bounding box";
        }
        Pout<< endl;

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

                Pout<< "    "
                    << setw(20) << pp.name()
                    << setw(9) << pp.size()
                    << setw(9) << pp.nPoints();


            primitivePatch::surfaceTopo pTyp = pp.surfaceType();

            if (pp.empty())
            {
                Pout<< setw(34) << "ok (empty)";
            }
            else if (pTyp == primitivePatch::MANIFOLD)
            {
                if (pp.checkPointManifold(true, &points))
                {
                    Pout<< setw(34) << "multiply connected (shared point)";
                }
                else
                {
                    Pout<< setw(34) << "ok (closed singly connected)";
                }

                // Add points on non-manifold edges to make set complete
                pp.checkTopology(false, &points);
            }
            else
            {
                pp.checkTopology(false, &points);

                if (pTyp == primitivePatch::OPEN)
                {
                    Pout<< setw(34) << "ok (non-closed singly connected)";
                }
                else
                {
                    Pout<< setw(34) << "multiply connected (shared edge)";
                }
            }

            if (allGeometry)
            {
                const pointField& pts = pp.points();
                const labelList& mp = pp.meshPoints();

                boundBox bb;   // zero-sized
                if (returnReduce(mp.size(), sumOp<label>()) > 0)
                {
                    bb.min() = pts[mp[0]];
                    bb.max() = pts[mp[0]];
                    for (label i = 1; i < mp.size(); i++)
                    {
                        bb.min() = min(bb.min(), pts[mp[i]]);
                        bb.max() = max(bb.max(), pts[mp[i]]);
                    }
                    reduce(bb.min(), minOp<vector>());
                    reduce(bb.max(), maxOp<vector>());
                }
                Pout<< ' ' << bb;
            }
            Pout<< endl;
        }

        if (points.size())
        {
            Pout<< "  <<Writing " << points.size()
                << " conflicting points to set "
                << points.name() << endl;

            points.write();
        }

        //Pout.setf(ios_base::right);
    }

    // Force creation of all addressing if requested.
    // Errors will be reported as required
    if (allTopology)
    {
        mesh.cells();
        mesh.faces();
        mesh.edges();
        mesh.points();
        mesh.faceOwner();
        mesh.faceNeighbour();
        mesh.cellCells();
        mesh.edgeCells();
        mesh.pointCells();
        mesh.edgeFaces();
        mesh.pointFaces();
        mesh.cellEdges();
        mesh.faceEdges();
        mesh.pointEdges();
    }

    return noFailedChecks;
}
