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

#include "meshSearch.H"
#include "polyMesh.H"
#include "indexedOctree.H"
#include "DynamicList.H"
#include "demandDrivenData.H"
#include "treeDataCell.H"
#include "treeDataFace.H"
#include "treeDataPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::meshSearch, 0);

Foam::scalar Foam::meshSearch::tol_ = 1E-3;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::meshSearch::findNearer
(
    const point& sample,
    const pointField& points,
    label& nearestI,
    scalar& nearestDistSqr
)
{
    bool nearer = false;

    forAll(points, pointI)
    {
        scalar distSqr = magSqr(points[pointI] - sample);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            nearestI = pointI;
            nearer = true;
        }
    }

    return nearer;
}


bool Foam::meshSearch::findNearer
(
    const point& sample,
    const pointField& points,
    const labelList& indices,
    label& nearestI,
    scalar& nearestDistSqr
)
{
    bool nearer = false;

    forAll(indices, i)
    {
        label pointI = indices[i];

        scalar distSqr = magSqr(points[pointI] - sample);

        if (distSqr < nearestDistSqr)
        {
            nearestDistSqr = distSqr;
            nearestI = pointI;
            nearer = true;
        }
    }

    return nearer;
}


// tree based searching
Foam::label Foam::meshSearch::findNearestCellTree(const point& location) const
{
    const indexedOctree<treeDataPoint>& tree = cellCentreTree();

    scalar span = tree.bb().mag();
    
    pointIndexHit info = tree.findNearest(location, Foam::sqr(span));

    if (!info.hit())
    {
        info = tree.findNearest(location, Foam::sqr(GREAT));
    }
    return info.index();
}


// linear searching
Foam::label Foam::meshSearch::findNearestCellLinear(const point& location) const
{
    const vectorField& centres = mesh_.cellCentres();

    label nearestIndex = 0;
    scalar minProximity = magSqr(centres[nearestIndex] - location);

    findNearer
    (
        location,
        centres,
        nearestIndex,
        minProximity
    );

    return nearestIndex;
}


// walking from seed
Foam::label Foam::meshSearch::findNearestCellWalk
(
    const point& location,
    const label seedCellI
) const
{
    if (seedCellI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findNearestCellWalk(const point&, const label)"
        )   << "illegal seedCell:" << seedCellI << exit(FatalError);
    }

    // Walk in direction of face that decreases distance

    label curCellI = seedCellI;
    scalar distanceSqr = magSqr(mesh_.cellCentres()[curCellI] - location);

    bool closer;

    do
    {
        // Try neighbours of curCellI
        closer = findNearer
        (
            location,
            mesh_.cellCentres(),
            mesh_.cellCells()[curCellI],
            curCellI,
            distanceSqr
        );
    } while (closer);

    return curCellI;
}


// tree based searching
Foam::label Foam::meshSearch::findNearestFaceTree(const point& location) const
{
    // Search nearest cell centre.
    const indexedOctree<treeDataPoint>& tree = cellCentreTree();

    scalar span = tree.bb().mag();

    // Search with decent span
    pointIndexHit info = tree.findNearest(location, Foam::sqr(span));

    if (!info.hit())
    {
        // Search with desparate span
        info = tree.findNearest(location, Foam::sqr(GREAT));
    }


    // Now check any of the faces of the nearest cell
    const vectorField& centres = mesh_.faceCentres();
    const cell& ownFaces = mesh_.cells()[info.index()];

    label nearestFaceI = ownFaces[0];
    scalar minProximity = magSqr(centres[nearestFaceI] - location);

    findNearer
    (
        location,
        centres,
        ownFaces,
        nearestFaceI,
        minProximity
    );

    return nearestFaceI;
}


// linear searching
Foam::label Foam::meshSearch::findNearestFaceLinear(const point& location) const
{
    const vectorField& centres = mesh_.faceCentres();

    label nearestFaceI = 0;
    scalar minProximity = magSqr(centres[nearestFaceI] - location);

    findNearer
    (
        location,
        centres,
        nearestFaceI,
        minProximity
    );

    return nearestFaceI;
}


// walking from seed
Foam::label Foam::meshSearch::findNearestFaceWalk
(
    const point& location,
    const label seedFaceI
) const
{
    if (seedFaceI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findNearestFaceWalk(const point&, const label)"
        )   << "illegal seedFace:" << seedFaceI << exit(FatalError);
    }

    const vectorField& centres = mesh_.faceCentres();


    // Walk in direction of face that decreases distance

    label curFaceI = seedFaceI;
    scalar distanceSqr = magSqr(centres[curFaceI] - location);

    while (true)
    {
        label betterFaceI = curFaceI;

        findNearer
        (
            location,
            centres,
            mesh_.cells()[mesh_.faceOwner()[curFaceI]],
            betterFaceI,
            distanceSqr
        );

        if (mesh_.isInternalFace(curFaceI))
        {
            findNearer
            (
                location,
                centres,
                mesh_.cells()[mesh_.faceNeighbour()[curFaceI]],
                betterFaceI,
                distanceSqr
            );
        }

        if (betterFaceI == curFaceI)
        {
            break;
        }

        curFaceI = betterFaceI;
    }

    return curFaceI;
}


Foam::label Foam::meshSearch::findCellLinear(const point& location) const
{
    bool cellFound = false;
    label n = 0;

    label cellI = -1;

    while ((!cellFound) && (n < mesh_.nCells()))
    {
        if (pointInCell(location, n))
        {
            cellFound = true;
            cellI = n;
        }
        else
        {
            n++;
        }
    }
    if (cellFound)
    {
        return cellI;
    }
    else
    {
        return -1;
    }
}


Foam::label Foam::meshSearch::findNearestBoundaryFaceWalk
(
    const point& location,
    const label seedFaceI
) const
{
    if (seedFaceI < 0)
    {
        FatalErrorIn
        (
            "meshSearch::findNearestBoundaryFaceWalk"
            "(const point&, const label)"
        )   << "illegal seedFace:" << seedFaceI << exit(FatalError);
    }

    // Start off from seedFaceI

    label curFaceI = seedFaceI;

    const face& f =  mesh_.faces()[curFaceI];

    scalar minDist = f.nearestPoint
    (
        location,
        mesh_.points()
    ).distance();

    bool closer;

    do
    {
        closer = false;

        // Search through all neighbouring boundary faces by going
        // across edges

        label lastFaceI = curFaceI;

        const labelList& myEdges = mesh_.faceEdges()[curFaceI];

        forAll (myEdges, myEdgeI)
        {
            const labelList& neighbours = mesh_.edgeFaces()[myEdges[myEdgeI]];

            // Check any face which uses edge, is boundary face and
            // is not curFaceI itself.

            forAll(neighbours, nI)
            {
                label faceI = neighbours[nI];

                if
                (
                    (faceI >= mesh_.nInternalFaces())
                 && (faceI != lastFaceI)
                )
                {
                    const face& f =  mesh_.faces()[faceI];

                    pointHit curHit = f.nearestPoint
                    (
                        location,
                        mesh_.points()
                    );

                    // If the face is closer, reset current face and distance
                    if (curHit.distance() < minDist)
                    {
                        minDist = curHit.distance();
                        curFaceI = faceI;
                        closer = true;  // a closer neighbour has been found
                    }
                }
            }
        }
    } while (closer);

    return curFaceI;
}


Foam::vector Foam::meshSearch::offset
(
    const point& bPoint,
    const label bFaceI,
    const vector& dir
) const
{
    // Get the neighbouring cell
    label ownerCellI = mesh_.faceOwner()[bFaceI];

    const point& c = mesh_.cellCentres()[ownerCellI];

    // Typical dimension: distance from point on face to cell centre
    scalar typDim = mag(c - bPoint);

    return tol_*typDim*dir;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshSearch::meshSearch(const polyMesh& mesh, const bool faceDecomp)
:
    mesh_(mesh),
    faceDecomp_(faceDecomp),
    cloud_(mesh_, IDLList<passiveParticle>()),
    boundaryTreePtr_(NULL),
    cellTreePtr_(NULL),
    cellCentreTreePtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshSearch::~meshSearch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataFace>& Foam::meshSearch::boundaryTree()
 const
{
    if (!boundaryTreePtr_)
    {
        //
        // Construct tree
        //

        // all boundary faces (not just walls)
        labelList bndFaces(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(bndFaces, i)
        {
            bndFaces[i] = mesh_.nInternalFaces() + i;
        }

        treeBoundBox overallBb(mesh_.points());
        Random rndGen(123456);
        overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        boundaryTreePtr_ = new indexedOctree<treeDataFace>
        (
            treeDataFace    // all information needed to search faces
            (
                false,                      // do not cache bb
                mesh_,
                bndFaces                    // boundary faces only
            ),
            overallBb,                      // overall search domain
            8,                              // maxLevel
            10,                             // leafsize
            3.0                             // duplicity
        );
    }

    return *boundaryTreePtr_;
}


const Foam::indexedOctree<Foam::treeDataCell>& Foam::meshSearch::cellTree()
 const
{
    if (!cellTreePtr_)
    {
        //
        // Construct tree
        //

        treeBoundBox overallBb(mesh_.points());
        Random rndGen(123456);
        overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        cellTreePtr_ = new indexedOctree<treeDataCell>
        (
            treeDataCell
            (
                false,  // not cache bb
                mesh_
            ),
            overallBb,  // overall search domain
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        );
    }

    return *cellTreePtr_;
    
}


const Foam::indexedOctree<Foam::treeDataPoint>&
 Foam::meshSearch::cellCentreTree() const
{
    if (!cellCentreTreePtr_)
    {
        //
        // Construct tree
        //

        treeBoundBox overallBb(mesh_.cellCentres());
        Random rndGen(123456);
        overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        cellCentreTreePtr_ = new indexedOctree<treeDataPoint>
        (
            treeDataPoint(mesh_.cellCentres()),
            overallBb,  // overall search domain
            10,         // max levels
            10.0,       // maximum ratio of cubes v.s. cells
            100.0       // max. duplicity; n/a since no bounding boxes.
        );
    }

    return *cellCentreTreePtr_;
}


// Is the point in the cell
// Works by checking if there is a face inbetween the point and the cell
// centre.
// Check for internal uses proper face decomposition or just average normal.
bool Foam::meshSearch::pointInCell(const point& p, label cellI) const
{
    if (faceDecomp_)
    {
        const point& ctr = mesh_.cellCentres()[cellI];

        vector dir(p - ctr);
        scalar magDir = mag(dir);

        // Check if any faces are hit by ray from cell centre to p.
        // If none -> p is in cell.
        const labelList& cFaces = mesh_.cells()[cellI];

        // Make sure half_ray does not pick up any faces on the wrong
        // side of the ray.
        scalar oldTol = intersection::setPlanarTol(0.0);

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            const face& f = mesh_.faces()[faceI];

            forAll(f, fp)
            {
                pointHit inter = f.ray
                (
                    ctr,
                    dir,
                    mesh_.points(),
                    intersection::HALF_RAY,
                    intersection::VECTOR
                );

                if (inter.hit())
                {
                    scalar dist = inter.distance();

                    if (dist < magDir)
                    {
                        // Valid hit. Hit face so point is not in cell.
                        intersection::setPlanarTol(oldTol);

                        return false;
                    }
                }
            }
        }

        intersection::setPlanarTol(oldTol);

        // No face inbetween point and cell centre so point is inside.
        return true;
    }
    else
    {
        const labelList& f = mesh_.cells()[cellI];
        const labelList& owner = mesh_.faceOwner();
        const vectorField& cf = mesh_.faceCentres();
        const vectorField& Sf = mesh_.faceAreas();

        forAll(f, facei)
        {
            label nFace = f[facei];
            vector proj = p - cf[nFace];
            vector normal = Sf[nFace];
            if (owner[nFace] == cellI)
            {
                if ((normal & proj) > 0)
                {
                    return false;
                }
            }
            else
            {
                if ((normal & proj) < 0)
                {
                    return false;
                }
            }
        }

        return true;
    }
}


Foam::label Foam::meshSearch::findNearestCell
(
    const point& location,
    const label seedCellI,
    const bool useTreeSearch
) const
{
    if (seedCellI == -1)
    {
        if (useTreeSearch)
        {
            return findNearestCellTree(location);
        }
        else
        {
            return findNearestCellLinear(location);
        }
    }
    else
    {
        return findNearestCellWalk(location, seedCellI);
    }
}


Foam::label Foam::meshSearch::findNearestFace
(
    const point& location,
    const label seedFaceI,
    const bool useTreeSearch
) const
{
    if (seedFaceI == -1)
    {
        if (useTreeSearch)
        {
            return findNearestFaceTree(location);
        }
        else
        {
            return findNearestFaceLinear(location);
        }
    }
    else
    {
        return findNearestFaceWalk(location, seedFaceI);
    }
}


Foam::label Foam::meshSearch::findCell
(
    const point& location,
    const label seedCellI,
    const bool useTreeSearch
) const
{
    // Find the nearest cell centre to this location
    label nearCellI = findNearestCell(location, seedCellI, useTreeSearch);

    if (debug)
    {
        Pout<< "findCell : nearCellI:" << nearCellI
            << " ctr:" << mesh_.cellCentres()[nearCellI]
            << endl;
    }

    if (pointInCell(location, nearCellI))
    {
        return nearCellI;
    }
    else
    {
        if (useTreeSearch)
        {
            // Start searching from cell centre of nearCell
            point curPoint = mesh_.cellCentres()[nearCellI];

            vector edgeVec = location - curPoint;
            edgeVec /= mag(edgeVec);

            do
            {
                // Walk neighbours (using tracking) to get closer
                passiveParticle tracker(cloud_, curPoint, nearCellI);

                if (debug)
                {
                    Pout<< "findCell : tracked from curPoint:" << curPoint
                        << " nearCellI:" << nearCellI;
                }

                tracker.track(location);

                if (debug)
                {
                    Pout<< " to " << tracker.position()
                        << " need:" << location
                        << " onB:" << tracker.onBoundary()
                        << " cell:" << tracker.cell()
                        << " face:" << tracker.face() << endl;
                }

                if (!tracker.onBoundary())
                {
                    // stopped not on boundary -> reached location
                    return tracker.cell();
                }

                // stopped on boundary face. Compare positions
                scalar typDim = sqrt(mag(mesh_.faceAreas()[tracker.face()]));

                if ((mag(tracker.position() - location)/typDim) < SMALL)
                {
                    return tracker.cell();
                }

                // tracking stopped at boundary. Find next boundary
                // intersection from current point (shifted out a little bit)
                // in the direction of the destination

                curPoint =
                    tracker.position()
                  + offset(tracker.position(), tracker.face(), edgeVec);

                if (debug)
                {
                    Pout<< "Searching for next boundary from curPoint:"
                        << curPoint
                        << " to location:" << location  << endl;
                }
                pointIndexHit curHit = intersection(curPoint, location);
                if (debug)
                {
                    Pout<< "Returned from line search with ";
                    curHit.write(Pout);
                    Pout<< endl;
                }

                if (!curHit.hit())
                {
                    return -1;
                }
                else
                {
                    nearCellI = mesh_.faceOwner()[curHit.index()];
                    curPoint =
                        curHit.hitPoint() 
                      + offset(curHit.hitPoint(), curHit.index(), edgeVec);
                }
            }
            while(true);
        }
        else
        {
             return findCellLinear(location);
        }
    }
    return -1;
}


Foam::label Foam::meshSearch::findNearestBoundaryFace
(
    const point& location,
    const label seedFaceI,
    const bool useTreeSearch
) const
{
    if (seedFaceI == -1)
    {
        if (useTreeSearch)
        {
            const indexedOctree<treeDataFace>& tree =  boundaryTree();

            scalar span = tree.bb().mag();

            pointIndexHit info = boundaryTree().findNearest
            (
                location,
                Foam::sqr(span)
            );

            if (!info.hit())
            {
                info = boundaryTree().findNearest
                (
                    location,
                    Foam::sqr(GREAT)
                );
            }

            return tree.shapes().faceLabels()[info.index()];
        }
        else
        {
            scalar minDist = GREAT;

            label minFaceI = -1;

            for
            (
                label faceI = mesh_.nInternalFaces();
                faceI < mesh_.nFaces();
                faceI++
            )
            {
                const face& f =  mesh_.faces()[faceI];

                pointHit curHit =
                    f.nearestPoint
                    (
                        location,
                        mesh_.points()
                    );

                if (curHit.distance() < minDist)
                {
                    minDist = curHit.distance();
                    minFaceI = faceI;
                }
            }
            return minFaceI;
        }
    }
    else
    {
        return findNearestBoundaryFaceWalk(location, seedFaceI);
    }
}


Foam::pointIndexHit Foam::meshSearch::intersection
(
    const point& pStart,
    const point& pEnd
) const
{
    pointIndexHit curHit = boundaryTree().findLine(pStart, pEnd);

    if (curHit.hit())
    {
        // Change index into octreeData into face label
        curHit.setIndex(boundaryTree().shapes().faceLabels()[curHit.index()]);
    }
    return curHit;
}


Foam::List<Foam::pointIndexHit> Foam::meshSearch::intersections
(
    const point& pStart,
    const point& pEnd
) const
{
    DynamicList<pointIndexHit> hits;

    vector edgeVec = pEnd - pStart;
    edgeVec /= mag(edgeVec);

    point pt = pStart;

    pointIndexHit bHit;
    do
    {
        bHit = intersection(pt, pEnd);

        if (bHit.hit())
        {
            hits.append(bHit);

            const vector& area = mesh_.faceAreas()[bHit.index()];

            scalar typDim = Foam::sqrt(mag(area));

            if ((mag(bHit.hitPoint() - pEnd)/typDim) < SMALL)
            {
                break;
            }

            // Restart from hitPoint shifted a little bit in the direction
            // of the destination

            pt =
                bHit.hitPoint()
              + offset(bHit.hitPoint(), bHit.index(), edgeVec);
        }

    } while(bHit.hit());


    hits.shrink();

    return hits;
}


bool Foam::meshSearch::isInside(const point& p) const
{
    return
        boundaryTree().getVolumeType(p)
     == indexedOctree<treeDataFace>::INSIDE;
}


// Delete all storage
void Foam::meshSearch::clearOut()
{
    deleteDemandDrivenData(boundaryTreePtr_);
    deleteDemandDrivenData(cellTreePtr_);
    deleteDemandDrivenData(cellCentreTreePtr_);
}


void Foam::meshSearch::correct()
{
    clearOut();
}


// ************************************************************************* //
