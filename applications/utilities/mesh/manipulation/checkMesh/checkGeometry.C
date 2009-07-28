#include "checkGeometry.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "EdgeMap.H"

Foam::label Foam::checkGeometry(const polyMesh& mesh, const bool allGeometry)
{
    label noFailedChecks = 0;

    Info<< "\nChecking geometry..." << endl;

    // Get a small relative length from the bounding box
    const boundBox& globalBb = mesh.bounds();

    Info<< "    Overall domain bounding box "
        << globalBb.min() << " " << globalBb.max() << endl;


    // Min length
    scalar minDistSqr = magSqr(1e-6 * globalBb.span());

    // Non-empty directions
    const Vector<label> validDirs = (mesh.geometricD() + Vector<label>::one)/2;
    Info<< "    Mesh (non-empty, non-wedge) directions " << validDirs << endl;

    const Vector<label> solDirs = (mesh.solutionD() + Vector<label>::one)/2;
    Info<< "    Mesh (non-empty) directions " << solDirs << endl;

    if (mesh.nGeometricD() < 3)
    {
        pointSet nonAlignedPoints(mesh, "nonAlignedEdges", mesh.nPoints()/100);

        if (mesh.checkEdgeAlignment(true, validDirs, &nonAlignedPoints))
        {
            noFailedChecks++;
            label nNonAligned = returnReduce
            (
                nonAlignedPoints.size(),
                sumOp<label>()
            );

            if (nNonAligned > 0)
            {
                Info<< "  <<Writing " << nNonAligned
                    << " points on non-aligned edges to set "
                    << nonAlignedPoints.name() << endl;
                nonAlignedPoints.write();
            }
        }
    }

    if (mesh.checkClosedBoundary(true)) noFailedChecks++;

    {
        cellSet cells(mesh, "nonClosedCells", mesh.nCells()/100+1);
        cellSet aspectCells(mesh, "highAspectRatioCells", mesh.nCells()/100+1);
        if (mesh.checkClosedCells(true, &cells, &aspectCells))
        {
            noFailedChecks++;

            label nNonClosed = returnReduce(cells.size(), sumOp<label>());

            if (nNonClosed > 0)
            {
                Info<< "  <<Writing " << nNonClosed
                    << " non closed cells to set " << cells.name() << endl;
                cells.write();
            }
        }

        label nHighAspect = returnReduce(aspectCells.size(), sumOp<label>());

        if (nHighAspect > 0)
        {
            Info<< "  <<Writing " << nHighAspect
                << " cells with high aspect ratio to set "
                << aspectCells.name() << endl;
            aspectCells.write();
        }
    }

    {
        faceSet faces(mesh, "zeroAreaFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceAreas(true, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " zero area faces to set " << faces.name() << endl;
                faces.write();
            }
        }
    }

    {
        cellSet cells(mesh, "zeroVolumeCells", mesh.nCells()/100 + 1);
        if (mesh.checkCellVolumes(true, &cells))
        {
            noFailedChecks++;

            label nCells = returnReduce(cells.size(), sumOp<label>());

            if (nCells > 0)
            {
                Info<< "  <<Writing " << nCells
                    << " zero volume cells to set " << cells.name() << endl;
                cells.write();
            }
        }
    }

    {
        faceSet faces(mesh, "nonOrthoFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceOrthogonality(true, &faces))
        {
            noFailedChecks++;
        }

        label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            Info<< "  <<Writing " << nFaces
                << " non-orthogonal faces to set " << faces.name() << endl;
            faces.write();
        }
    }


    {
        faceSet faces(mesh, "wrongOrientedFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFacePyramids(true, -SMALL, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " faces with incorrect orientation to set "
                    << faces.name() << endl;
                faces.write();
            }
        }
    }

    {
        faceSet faces(mesh, "skewFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceSkewness(true, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " skew faces to set " << faces.name() << endl;
                faces.write();
            }
        }
    }

    if (allGeometry)
    {
        // Note use of nPoints since don't want edge construction.
        pointSet points(mesh, "shortEdges", mesh.nPoints()/1000 + 1);
        if (mesh.checkEdgeLength(true, minDistSqr, &points))
        {
            //noFailedChecks++;

            label nPoints = returnReduce(points.size(), sumOp<label>());

            if (nPoints > 0)
            {
                Info<< "  <<Writing " << nPoints
                    << " points on short edges to set " << points.name()
                    << endl;
                points.write();
            }
        }

        label nEdgeClose = returnReduce(points.size(), sumOp<label>());

        if (mesh.checkPointNearness(false, minDistSqr, &points))
        {
            //noFailedChecks++;

            label nPoints = returnReduce(points.size(), sumOp<label>());

            if (nPoints > nEdgeClose)
            {
                pointSet nearPoints(mesh, "nearPoints", points);
                Info<< "  <<Writing " << nPoints
                    << " near (closer than " << Foam::sqrt(minDistSqr)
                    << " apart) points to set " << nearPoints.name() << endl;
                nearPoints.write();
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "concaveFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceAngles(true, 10, &faces))
        {
            //noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " faces with concave angles to set " << faces.name()
                    << endl;
                faces.write();
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "warpedFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceFlatness(true, 0.8, &faces))
        {
            //noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " warped faces to set " << faces.name() << endl;
                faces.write();
            }
        }
    }

    if (allGeometry)
    {
        cellSet cells(mesh, "underdeterminedCells", mesh.nCells()/100);
        if (mesh.checkCellDeterminant(true, &cells))
        {
            noFailedChecks++;

            label nCells = returnReduce(cells.size(), sumOp<label>());

            Info<< "  <<Writing " << nCells
                << " under-determined cells to set " << cells.name() << endl;
            cells.write();
        }
    }


    return noFailedChecks;
}
