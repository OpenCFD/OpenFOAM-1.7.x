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

#include "cyclicPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "globalPointPatch.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicPointPatch, 0);

addToRunTimeSelectionTable
(
    facePointPatch,
    cyclicPointPatch,
    polyPatch
);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicPointPatch::initGeometry()
{
    transformPairs_.setSize(0);
}


void Foam::cyclicPointPatch::calcGeometry()
{
    const edgeList& cp = cyclicPolyPatch_.coupledPoints();
    const labelList& mp = cyclicPolyPatch_.meshPoints();

    // If there are no global points create a 1->1 map
    if (!boundaryMesh().mesh().globalData().nGlobalPoints())
    {
        nonGlobalPatchPoints_.setSize(mp.size());
        forAll(nonGlobalPatchPoints_, i)
        {
            nonGlobalPatchPoints_[i] = i;
        }

        meshPoints_ = cyclicPolyPatch_.meshPoints();
        transformPairs_ = cp;
    }
    else
    {
        // Get reference to shared points
        const labelList& sharedPoints = 
            boundaryMesh().globalPatch().meshPoints();

        nonGlobalPatchPoints_.setSize(mp.size());
        meshPoints_.setSize(mp.size());

        labelList pointMap(mp.size(), -1);

        label noFiltPoints = 0;

        forAll (mp, pointI)
        {
            label curP = mp[pointI];

            bool found = false;

            forAll (sharedPoints, sharedI)
            {
                if (sharedPoints[sharedI] == curP)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                pointMap[pointI] = noFiltPoints;
                nonGlobalPatchPoints_[noFiltPoints] = pointI;
                meshPoints_[noFiltPoints] = curP;
                noFiltPoints++;
            }
        }

        nonGlobalPatchPoints_.setSize(noFiltPoints);
        meshPoints_.setSize(noFiltPoints);


        transformPairs_.setSize(cp.size());

        label noFiltPointPairs = 0;

        forAll(cp, i)
        {
            if (pointMap[cp[i][0]] != -1 && pointMap[cp[i][1]] != -1)
            {
                transformPairs_[noFiltPointPairs][0] = pointMap[cp[i][0]];
                transformPairs_[noFiltPointPairs][1] = pointMap[cp[i][1]];
                noFiltPointPairs++;
            }
            else if (pointMap[cp[i][0]] == -1 && pointMap[cp[i][1]] != -1)
            {
                FatalErrorIn("cyclicPointPatch::calcGeometry() const")
                    << "Point " << cp[i][0] << "of point-pair " << i
                    << " is a global point but the other point "
                    << cp[i][1] << " is not"
                    << exit(FatalError);
            }
            else if (pointMap[cp[i][0]] != -1 && pointMap[cp[i][1]] == -1)
            {
                FatalErrorIn("cyclicPointPatch::calcGeometry() const")
                    << "Point " << cp[i][1] << "of point-pair " << i
                    << " is a global point but the other point "
                    << cp[i][0] << " is not"
                    << exit(FatalError);
            }
        }

        transformPairs_.setSize(noFiltPointPairs);
    }
}


void cyclicPointPatch::initMovePoints(const pointField&)
{}


void cyclicPointPatch::movePoints(const pointField&)
{}


void cyclicPointPatch::initUpdateMesh()
{
    facePointPatch::initUpdateMesh();
    cyclicPointPatch::initGeometry();
}


void cyclicPointPatch::updateMesh()
{
    facePointPatch::updateMesh();
    cyclicPointPatch::calcGeometry();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cyclicPointPatch::cyclicPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    cyclicPolyPatch_(refCast<const cyclicPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cyclicPointPatch::~cyclicPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const edgeList& cyclicPointPatch::transformPairs() const
{
    return transformPairs_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
