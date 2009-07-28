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

#include "domainDecomposition.H"
#include "Time.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "fvMesh.H"
#include "OSspecific.H"
#include "Map.H"
#include "globalMeshData.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void domainDecomposition::mark
(
    const labelList& zoneElems,
    const label zoneI,
    labelList& elementToZone
)
{
    forAll(zoneElems, i)
    {
        label pointi = zoneElems[i];

        if (elementToZone[pointi] == -1)
        {
            // First occurrence
            elementToZone[pointi] = zoneI;
        }
        else if (elementToZone[pointi] >= 0)
        {
            // Multiple zones
            elementToZone[pointi] = -2;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
domainDecomposition::domainDecomposition(const IOobject& io)
:
    fvMesh(io),
    decompositionDict_
    (
        IOobject
        (
            "decomposeParDict",
            time().system(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nProcs_(readInt(decompositionDict_.lookup("numberOfSubdomains"))),
    distributed_(false),
    cellToProc_(nCells()),
    procPointAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    procCellAddressing_(nProcs_),
    procBoundaryAddressing_(nProcs_),
    procPatchSize_(nProcs_),
    procPatchStartIndex_(nProcs_),
    procNeighbourProcessors_(nProcs_),
    procProcessorPatchSize_(nProcs_),
    procProcessorPatchStartIndex_(nProcs_),
    globallySharedPoints_(0),
    cyclicParallel_(false)
{
    if (decompositionDict_.found("distributed"))
    {
        Switch distributed(decompositionDict_.lookup("distributed"));
        distributed_ = distributed;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

domainDecomposition::~domainDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool domainDecomposition::writeDecomposition()
{
    Info<< "\nConstructing processor meshes" << endl;

    // Make a lookup map for globally shared points
    Map<label> sharedPointLookup(2*globallySharedPoints_.size());

    forAll (globallySharedPoints_, pointi)
    {
        sharedPointLookup.insert(globallySharedPoints_[pointi], pointi);
    }


    // Mark point/faces/cells that are in zones.
    // -1   : not in zone
    // -2   : in multiple zones
    // >= 0 : in single given zone
    // This will give direct lookup of elements that are in a single zone
    // and we'll only have to revert back to searching through all zones
    // for the duplicate elements

    // Point zones
    labelList pointToZone(points().size(), -1);

    forAll(pointZones(), zoneI)
    {
        mark(pointZones()[zoneI], zoneI, pointToZone);
    }

    // Face zones
    labelList faceToZone(faces().size(), -1);

    forAll(faceZones(), zoneI)
    {
        mark(faceZones()[zoneI], zoneI, faceToZone);
    }

    // Cell zones
    labelList cellToZone(nCells(), -1);

    forAll(cellZones(), zoneI)
    {
        mark(cellZones()[zoneI], zoneI, cellToZone);
    }


    label totProcFaces = 0;
    label maxProcPatches = 0;
    label maxProcFaces = 0;


    // Write out the meshes
    for (label procI = 0; procI < nProcs_; procI++)
    {
        // Create processor points
        const labelList& curPointLabels = procPointAddressing_[procI];

        const pointField& meshPoints = points();

        labelList pointLookup(nPoints(), -1);

        pointField procPoints(curPointLabels.size());

        forAll (curPointLabels, pointi)
        {
            procPoints[pointi] = meshPoints[curPointLabels[pointi]];

            pointLookup[curPointLabels[pointi]] = pointi;
        }

        // Create processor faces
        const labelList& curFaceLabels = procFaceAddressing_[procI];

        const faceList& meshFaces = faces();

        labelList faceLookup(nFaces(), -1);

        faceList procFaces(curFaceLabels.size());

        forAll (curFaceLabels, facei)
        {
            // Mark the original face as used
            // Remember to decrement the index by one (turning index)
            //
            label curF = mag(curFaceLabels[facei]) - 1;

            faceLookup[curF] = facei;

            // get the original face
            labelList origFaceLabels;

            if (curFaceLabels[facei] >= 0)
            {
                // face not turned
                origFaceLabels = meshFaces[curF];
            }
            else
            {
                origFaceLabels = meshFaces[curF].reverseFace();
            }

            // translate face labels into local point list
            face& procFaceLabels = procFaces[facei];

            procFaceLabels.setSize(origFaceLabels.size());

            forAll (origFaceLabels, pointi)
            {
                procFaceLabels[pointi] = pointLookup[origFaceLabels[pointi]];
            }
        }

        // Create processor cells
        const labelList& curCellLabels = procCellAddressing_[procI];

        const cellList& meshCells = cells();

        cellList procCells(curCellLabels.size());

        forAll (curCellLabels, celli)
        {
            const labelList& origCellLabels = meshCells[curCellLabels[celli]];

            cell& curCell = procCells[celli];

            curCell.setSize(origCellLabels.size());

            forAll (origCellLabels, cellFaceI)
            {
                curCell[cellFaceI] = faceLookup[origCellLabels[cellFaceI]];
            }
        }

        // Create processor mesh without a boundary

        fileName processorCasePath
        (
            time().caseName()/fileName(word("processor") + Foam::name(procI))
        );

        // make the processor directory
        mkDir(time().rootPath()/processorCasePath);

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            processorCasePath,
            "system",
            "constant"
        );

        // create the mesh
        polyMesh procMesh
        (
            IOobject
            (
                this->polyMesh::name(),  // region name of undecomposed mesh
                pointsInstance(),
                processorDb
            ),
            xferMove(procPoints),
            xferMove(procFaces),
            xferMove(procCells)
        );

        // Create processor boundary patches
        const labelList& curBoundaryAddressing = procBoundaryAddressing_[procI];

        const labelList& curPatchSizes = procPatchSize_[procI];

        const labelList& curPatchStarts = procPatchStartIndex_[procI];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        const labelList& curProcessorPatchStarts =
            procProcessorPatchStartIndex_[procI];

        const polyPatchList& meshPatches = boundaryMesh();

        List<polyPatch*> procPatches
        (
            curPatchSizes.size()
          + curProcessorPatchSizes.size(),
            reinterpret_cast<polyPatch*>(0)
        );

        label nPatches = 0;

        forAll (curPatchSizes, patchi)
        {
            procPatches[nPatches] =
                meshPatches[curBoundaryAddressing[patchi]].clone
                (
                    procMesh.boundaryMesh(),
                    nPatches,
                    curPatchSizes[patchi],
                    curPatchStarts[patchi]
                ).ptr();

            nPatches++;
        }

        forAll (curProcessorPatchSizes, procPatchI)
        {
            procPatches[nPatches] =
                new processorPolyPatch
                (
                    word("procBoundary") + Foam::name(procI)
                  + word("to")
                  + Foam::name(curNeighbourProcessors[procPatchI]),
                    curProcessorPatchSizes[procPatchI],
                    curProcessorPatchStarts[procPatchI],
                    nPatches,
                    procMesh.boundaryMesh(),
                    procI,
                    curNeighbourProcessors[procPatchI]
            );

            nPatches++;
        }

        // Add boundary patches
        procMesh.addPatches(procPatches);

        // Create and add zones

        // Point zones
        {
            const pointZoneMesh& pz = pointZones();

            // Go through all the zoned points and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label> > zonePoints(pz.size());

            // Estimate size
            forAll(zonePoints, zoneI)
            {
                zonePoints[zoneI].setCapacity(pz[zoneI].size() / nProcs_);
            }

            // Use the pointToZone map to find out the single zone (if any),
            // use slow search only for shared points.
            forAll (curPointLabels, pointi)
            {
                label curPoint = curPointLabels[pointi];

                label zoneI = pointToZone[curPoint];

                if (zoneI >= 0)
                {
                    // Single zone.
                    zonePoints[zoneI].append(pointi);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(pz, zoneI)
                    {
                        label index = pz[zoneI].whichPoint(curPoint);

                        if (index != -1)
                        {
                            zonePoints[zoneI].append(pointi);
                        }
                    }
                }
            }

            procMesh.pointZones().clearAddressing();
            procMesh.pointZones().setSize(zonePoints.size());
            forAll(zonePoints, zoneI)
            {
                procMesh.pointZones().set
                (
                    zoneI,
                    pz[zoneI].clone
                    (
                        procMesh.pointZones(),
                        zoneI,
                        zonePoints[zoneI].shrink()
                    )
                );
            }

            if (pz.size())
            {
                // Force writing on all processors
                procMesh.pointZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Face zones
        {
            const faceZoneMesh& fz = faceZones();

            // Go through all the zoned face and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label> > zoneFaces(fz.size());
            List<DynamicList<bool> > zoneFaceFlips(fz.size());

            // Estimate size
            forAll(zoneFaces, zoneI)
            {
                label procSize = fz[zoneI].size() / nProcs_;

                zoneFaces[zoneI].setCapacity(procSize);
                zoneFaceFlips[zoneI].setCapacity(procSize);
            }

            // Go through all the zoned faces and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            forAll (curFaceLabels, facei)
            {
                // Remember to decrement the index by one (turning index)
                //
                label curF = mag(curFaceLabels[facei]) - 1;

                label zoneI = faceToZone[curF];

                if (zoneI >= 0)
                {
                    // Single zone. Add the face
                    zoneFaces[zoneI].append(facei);

                    label index = fz[zoneI].whichFace(curF);

                    bool flip = fz[zoneI].flipMap()[index];

                    if (curFaceLabels[facei] < 0)
                    {
                        flip = !flip;
                    }

                    zoneFaceFlips[zoneI].append(flip);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(fz, zoneI)
                    {
                        label index = fz[zoneI].whichFace(curF);

                        if (index != -1)
                        {
                            zoneFaces[zoneI].append(facei);

                            bool flip = fz[zoneI].flipMap()[index];

                            if (curFaceLabels[facei] < 0)
                            {
                                flip = !flip;
                            }

                            zoneFaceFlips[zoneI].append(flip);
                        }
                    }
                }
            }

            procMesh.faceZones().clearAddressing();
            procMesh.faceZones().setSize(zoneFaces.size());
            forAll(zoneFaces, zoneI)
            {
                procMesh.faceZones().set
                (
                    zoneI,
                    fz[zoneI].clone
                    (
                        zoneFaces[zoneI].shrink(),          // addressing
                        zoneFaceFlips[zoneI].shrink(),      // flipmap
                        zoneI,
                        procMesh.faceZones()
                    )
                );
            }

            if (fz.size())
            {
                // Force writing on all processors
                procMesh.faceZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Cell zones
        {
            const cellZoneMesh& cz = cellZones();

            // Go through all the zoned cells and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label> > zoneCells(cz.size());

            // Estimate size
            forAll(zoneCells, zoneI)
            {
                zoneCells[zoneI].setCapacity(cz[zoneI].size() / nProcs_);
            }

            forAll (curCellLabels, celli)
            {
                label curCellI = curCellLabels[celli];

                label zoneI = cellToZone[curCellI];

                if (zoneI >= 0)
                {
                    // Single zone.
                    zoneCells[zoneI].append(celli);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(cz, zoneI)
                    {
                        label index = cz[zoneI].whichCell(curCellI);

                        if (index != -1)
                        {
                            zoneCells[zoneI].append(celli);
                        }
                    }
                }
            }

            procMesh.cellZones().clearAddressing();
            procMesh.cellZones().setSize(zoneCells.size());
            forAll(zoneCells, zoneI)
            {
                procMesh.cellZones().set
                (
                    zoneI,
                    cz[zoneI].clone
                    (
                        zoneCells[zoneI].shrink(),
                        zoneI,
                        procMesh.cellZones()
                    )
                );
            }

            if (cz.size())
            {
                // Force writing on all processors
                procMesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(10);

        procMesh.write();

        Info<< endl
            << "Processor " << procI << nl
            << "    Number of cells = " << procMesh.nCells()
            << endl;

        label nBoundaryFaces = 0;
        label nProcPatches = 0;
        label nProcFaces = 0;

        forAll (procMesh.boundaryMesh(), patchi)
        {
            if
            (
                procMesh.boundaryMesh()[patchi].type()
             == processorPolyPatch::typeName
            )
            {
                const processorPolyPatch& ppp =
                refCast<const processorPolyPatch>
                (
                    procMesh.boundaryMesh()[patchi]
                );

                Info<< "    Number of faces shared with processor "
                    << ppp.neighbProcNo() << " = " << ppp.size() << endl;

                nProcPatches++;
                nProcFaces += ppp.size();
            }
            else
            {
                nBoundaryFaces += procMesh.boundaryMesh()[patchi].size();
            }
        }

        Info<< "    Number of processor patches = " << nProcPatches << nl
            << "    Number of processor faces = " << nProcFaces << nl
            << "    Number of boundary faces = " << nBoundaryFaces << endl;

        totProcFaces += nProcFaces;
        maxProcPatches = max(maxProcPatches, nProcPatches);
        maxProcFaces = max(maxProcFaces, nProcFaces);

        // create and write the addressing information
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procPointAddressing_[procI]
        );
        pointProcAddressing.write();

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procFaceAddressing_[procI]
        );
        faceProcAddressing.write();

        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procCellAddressing_[procI]
        );
        cellProcAddressing.write();

        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procBoundaryAddressing_[procI]
        );
        boundaryProcAddressing.write();
    }

    Info<< nl
        << "Number of processor faces = " << totProcFaces/2 << nl
        << "Max number of processor patches = " << maxProcPatches << nl
        << "Max number of faces between processors = " << maxProcFaces
        << endl;

    return true;
}


// ************************************************************************* //
