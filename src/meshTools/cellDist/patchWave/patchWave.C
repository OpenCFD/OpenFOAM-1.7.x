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

\*---------------------------------------------------------------------------*/

#include "patchWave.H"
#include "polyMesh.H"
#include "wallPoint.H"
#include "MeshWave.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchWave::setChangedFaces
(
    const labelHashSet& patchIDs,
    labelList& changedFaces,
    List<wallPoint>& faceDist
) const
{
    const polyMesh& mesh = cellDistFuncs::mesh();

    label nChangedFaces = 0;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchIDs.found(patchI))
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];

            forAll(patch.faceCentres(), patchFaceI)
            {
                label meshFaceI = patch.start() + patchFaceI;

                changedFaces[nChangedFaces] = meshFaceI;

                faceDist[nChangedFaces] =
                    wallPoint
                    (
                        patch.faceCentres()[patchFaceI],
                        0.0
                    );

                nChangedFaces++;
            }
        }
    }
}


Foam::label Foam::patchWave::getValues(const MeshWave<wallPoint>& waveInfo)
{
    const List<wallPoint>& cellInfo = waveInfo.allCellInfo();
    const List<wallPoint>& faceInfo = waveInfo.allFaceInfo();

    label nIllegal = 0;

    // Copy cell values
    distance_.setSize(cellInfo.size());

    forAll(cellInfo, cellI)
    {
        scalar dist = cellInfo[cellI].distSqr();

        if (cellInfo[cellI].valid())
        {
            distance_[cellI] = Foam::sqrt(dist);
        }
        else
        {
            distance_[cellI] = dist;

            nIllegal++;
        }
    }

    // Copy boundary values
    forAll(patchDistance_, patchI)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchI];

        // Allocate storage for patchDistance
        scalarField* patchDistPtr = new scalarField(patch.size());

        patchDistance_.set(patchI, patchDistPtr);

        scalarField& patchField = *patchDistPtr;

        forAll(patchField, patchFaceI)
        {
            label meshFaceI = patch.start() + patchFaceI;

            scalar dist = faceInfo[meshFaceI].distSqr();

            if (faceInfo[meshFaceI].valid())
            {
                // Adding SMALL to avoid problems with /0 in the turbulence
                // models
                patchField[patchFaceI] = Foam::sqrt(dist) + SMALL;
            }
            else
            {
                patchField[patchFaceI] = dist;

                nIllegal++;
            }
        }
    }
    return nIllegal;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::patchWave::patchWave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    const bool correctWalls
)
:
    cellDistFuncs(mesh),
    patchIDs_(patchIDs),
    correctWalls_(correctWalls),
    nUnset_(0),
    distance_(mesh.nCells()),
    patchDistance_(mesh.boundaryMesh().size())
{
    patchWave::correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchWave::~patchWave()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Correct for mesh geom/topo changes. Might be more intelligent in the
// future (if only small topology change)
void Foam::patchWave::correct()
{
    //
    // Set initial changed faces: set wallPoint for wall faces to wall centre
    //

    // Count walls
    label nWalls = sumPatchSize(patchIDs_);

    List<wallPoint> faceDist(nWalls);
    labelList changedFaces(nWalls);

    // Set to faceDist information to facecentre on walls.
    setChangedFaces(patchIDs_, changedFaces, faceDist);


    //
    // Do calculate wall distance by 'growing' from faces.
    //

    MeshWave<wallPoint> waveInfo
    (
        mesh(),
        changedFaces,
        faceDist,
        mesh().globalData().nTotalCells()   // max iterations
    );


    //
    // Copy distance into return field
    //

    nUnset_ = getValues(waveInfo);

    //
    // Correct wall cells for true distance
    //

    if (correctWalls_)
    {
        Map<label> nearestFace(2 * nWalls);

        correctBoundaryFaceCells
        (
            patchIDs_,
            distance_,
            nearestFace
        );

        correctBoundaryPointCells
        (
            patchIDs_,
            distance_,
            nearestFace
        );
    }
}


// ************************************************************************* //
