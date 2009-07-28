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

#include "processorPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "globalPointPatch.H"
#include "faceList.H"
#include "primitiveFacePatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorPointPatch, 0);

addToRunTimeSelectionTable
(
    facePointPatch,
    processorPointPatch,
    polyPatch
);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::processorPointPatch::initGeometry()
{
    // Algorithm:
    // Depending on whether the patch is a master or a slave, get the primitive
    // patch points and filter away the points from the global patch.

    if (isMaster())
    {
        meshPoints_ = procPolyPatch_.meshPoints();
    }
    else
    {
        // Slave side. Create the reversed patch and pick up its points
        // so that the order is correct
        const polyPatch& pp = patch();

        faceList masterFaces(pp.size());

        forAll (pp, faceI)
        {
            masterFaces[faceI] = pp[faceI].reverseFace();
        }

        meshPoints_ = primitiveFacePatch
        (
            masterFaces,
            pp.points()
        ).meshPoints();
    }

    if (Pstream::parRun())
    {
        initPatchPatchPoints();
    }
}


void Foam::processorPointPatch::calcGeometry()
{
    if (Pstream::parRun())
    {
        calcPatchPatchPoints();
    }

    // If it is not runing parallel or there are no global points
    // create a 1->1 map
    if
    (
        !Pstream::parRun()
     || !boundaryMesh().mesh().globalData().nGlobalPoints()
    )
    {
        nonGlobalPatchPoints_.setSize(meshPoints_.size());
        forAll(nonGlobalPatchPoints_, i)
        {
            nonGlobalPatchPoints_[i] = i;
        }
    }
    else
    {
        // Get reference to shared points
        const labelList& sharedPoints =
            boundaryMesh().globalPatch().meshPoints();

        nonGlobalPatchPoints_.setSize(meshPoints_.size());

        label noFiltPoints = 0;

        forAll (meshPoints_, pointI)
        {
            label curP = meshPoints_[pointI];

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
                nonGlobalPatchPoints_[noFiltPoints] = pointI;
                meshPoints_[noFiltPoints] = curP;
                noFiltPoints++;
            }
        }

        nonGlobalPatchPoints_.setSize(noFiltPoints);
        meshPoints_.setSize(noFiltPoints);
    }
}


void processorPointPatch::initPatchPatchPoints()
{
    if (debug)
    {
        Info<< "processorPointPatch::calcPatchPatchPoints() : "
            << "constructing patch-patch points"
            << endl;
    }

    const polyBoundaryMesh& bm = boundaryMesh().mesh()().boundaryMesh();

    // Get the mesh points for this patch corresponding to the faces
    const labelList& ppmp = meshPoints();

    // Create a HashSet of the point labels for this patch
    Map<label> patchPointSet(2*ppmp.size());

    forAll (ppmp, ppi)
    {
        patchPointSet.insert(ppmp[ppi], ppi);
    }


    // Create the lists of patch-patch points
    labelListList patchPatchPoints(bm.size());

    // Create the lists of patch-patch point normals
    List<List<vector> > patchPatchPointNormals(bm.size());

    // Loop over all patches looking for other patches that share points
    forAll(bm, patchi)
    {
        if
        (
            patchi != index()                 // Ignore self-self
         && !isA<emptyPolyPatch>(bm[patchi])  // Ignore empty
         && !bm[patchi].coupled()             // Ignore other couples
        )
        {
            // Get the meshPoints for the other patch
            const labelList& meshPoints = bm[patchi].meshPoints();

            // Get the normals for the other patch
            const vectorField& normals = bm[patchi].pointNormals();

            label pppi = 0;
            forAll(meshPoints, pointi)
            {
                label ppp = meshPoints[pointi];

                // Check to see if the point of the other patch is shared with
                // this patch
                Map<label>::iterator iter = patchPointSet.find(ppp);

                if (iter != patchPointSet.end())
                {
                    // If it is shared initialise the patchPatchPoints for this
                    // patch
                    if (!patchPatchPoints[patchi].size())
                    {
                        patchPatchPoints[patchi].setSize(ppmp.size());
                        patchPatchPointNormals[patchi].setSize(ppmp.size());
                    }

                    // and add the entry
                    patchPatchPoints[patchi][pppi] = iter();
                    patchPatchPointNormals[patchi][pppi] = normals[pointi];
                    pppi++;
                }
            }

            // Resise the list of shared points and normals for the patch
            // being considerd
            patchPatchPoints[patchi].setSize(pppi);
            patchPatchPointNormals[patchi].setSize(pppi);
        }
    }

    // Send the patchPatchPoints to the neighbouring processor

    OPstream toNeighbProc
    (
        Pstream::blocking,
        neighbProcNo()
    );

    toNeighbProc
        << ppmp.size()              // number of points for checking
        << patchPatchPoints
        << patchPatchPointNormals;

    if (debug)
    {
        Info<< "processorPointPatch::calcPatchPatchPoints() : "
            << "constructed patch-patch points"
            << endl;
    }
}


void Foam::processorPointPatch::calcPatchPatchPoints()
{
    // Get the patchPatchPoints from the neighbouring processor
    IPstream fromNeighbProc
    (
        Pstream::blocking,
        neighbProcNo()
    );

    label nbrNPoints(readLabel(fromNeighbProc));
    labelListList patchPatchPoints(fromNeighbProc);
    List<List<vector> > patchPatchPointNormals(fromNeighbProc);

    pointBoundaryMesh& pbm = const_cast<pointBoundaryMesh&>(boundaryMesh());
    const labelList& ppmp = meshPoints();

    // Simple check for the very rare situation when not the same number
    // of points on both sides. This can happen with decomposed cyclics.
    // If on one side the cyclic shares a point with proc faces coming from
    // internal faces it will have a different number of points from
    // the situation where the cyclic and the 'normal' proc faces are fully
    // separate.
    if (nbrNPoints != ppmp.size())
    {
        WarningIn("processorPointPatch::calcPatchPatchPoints()")
            << "Processor patch " << name()
            << " has " << ppmp.size() << " points; coupled patch has "
            << nbrNPoints << " points." << endl
            << "   (usually due to decomposed cyclics)."
            << " This might give problems" << endl
            << "    when using point fields (interpolation, mesh motion)."
            << endl;
    }



    // Loop over the patches looking for other patches that share points
    forAll(patchPatchPoints, patchi)
    {
        const labelList& patchPoints = patchPatchPoints[patchi];
        const List<vector>& patchPointNormals = patchPatchPointNormals[patchi];

        // If there are potentially shared points for the patch being considered
        if (patchPoints.size())
        {
            // Get the current meshPoints list for the patch
            facePointPatch& fpp = refCast<facePointPatch>(pbm[patchi]);
            const labelList& fmp = fpp.meshPoints();
            labelList& mp = fpp.meshPoints_;

            const vectorField& fnormals = fpp.pointNormals();
            vectorField& normals = fpp.pointNormals_;

            // Create a HashSet of the point labels for the patch
            Map<label> patchPointSet(2*fmp.size());

            forAll (fmp, ppi)
            {
                patchPointSet.insert(fmp[ppi], ppi);
            }

            label nPoints = mp.size();
            label lpi = 0;
            bool resized = false;

            // For each potentially shared point...
            forAll(patchPoints, ppi)
            {
                // Check if it is not already in the patch,
                // i.e. not part of a face of the patch
                if (!patchPointSet.found(ppmp[patchPoints[ppi]]))
                {
                    // If it isn't already in the patch check if the local
                    // meshPoints is already set and if not initialise the
                    // meshPoints_ and pointNormals_
                    if (!resized)
                    {
                        if (!mp.size() && fmp.size())
                        {
                            mp = fmp;
                            normals = fnormals;

                            nPoints = mp.size();
                        }

                        mp.setSize(nPoints + patchPoints.size());
                        loneMeshPoints_.setSize(patchPoints.size());
                        normals.setSize(nPoints + patchPoints.size());
                        resized = true;
                    }

                    // Add the new point to the patch
                    mp[nPoints] = ppmp[patchPoints[ppi]];
                    loneMeshPoints_[lpi++] = ppmp[patchPoints[ppi]];
                    normals[nPoints++] = patchPointNormals[ppi];
                }
            }

            // If the lists have been resized points have been added.
            // Shrink the lists to the current size.
            if (resized)
            {
                mp.setSize(nPoints);
                loneMeshPoints_.setSize(lpi);
                normals.setSize(nPoints);
            }
        }
    }
}


void processorPointPatch::initMovePoints(const pointField&)
{}


void processorPointPatch::movePoints(const pointField&)
{}


void processorPointPatch::initUpdateMesh()
{
    facePointPatch::initUpdateMesh();
    processorPointPatch::initGeometry();
}


void processorPointPatch::updateMesh()
{
    facePointPatch::updateMesh();
    processorPointPatch::calcGeometry();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

processorPointPatch::processorPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    procPolyPatch_(refCast<const processorPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorPointPatch::~processorPointPatch()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
