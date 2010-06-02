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

Description
    Utility to refine cells next to patches.

    Takes a patchName and number of layers to refine. Works out cells within
    these layers and refines those in the wall-normal direction.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "cellCuts.H"
#include "cellSet.H"
#include "meshCutter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::noParallel();
    Foam::argList::validArgs.append("patchName");
    Foam::argList::validArgs.append("edgeWeight");
    Foam::argList::validOptions.insert("useSet", "cellSet");
    Foam::argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    word patchName(args.additionalArgs()[0]);

    scalar weight(readScalar(IStringStream(args.additionalArgs()[1])()));
    bool overwrite = args.optionFound("overwrite");


    label patchID = mesh.boundaryMesh().findPatchID(patchName);

    if (patchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Cannot find patch " << patchName << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }
    const polyPatch& pp = mesh.boundaryMesh()[patchID];


    // Cells cut

    labelHashSet cutCells(4*pp.size());

    const labelList& meshPoints = pp.meshPoints();

    forAll(meshPoints, pointI)
    {
        label meshPointI = meshPoints[pointI];

        const labelList& pCells = mesh.pointCells()[meshPointI];

        forAll(pCells, pCellI)
        {
            cutCells.insert(pCells[pCellI]);
        }
    }

    Info<< "Selected " << cutCells.size()
        << " cells connected to patch " << pp.name() << endl << endl;

    //
    // List of cells to refine
    //

    bool useSet = args.optionFound("useSet");

    if (useSet)
    {
        word setName(args.option("useSet"));

        Info<< "Subsetting cells to cut based on cellSet" << setName << endl
            << endl;

        cellSet cells(mesh, setName);

        Info<< "Read " << cells.size() << " cells from cellSet "
            << cells.instance()/cells.local()/cells.name()
            << endl << endl;

        for
        (
            cellSet::const_iterator iter = cells.begin();
            iter != cells.end();
            ++iter
        )
        {
            cutCells.erase(iter.key());
        }
        Info<< "Removed from cells to cut all the ones not in set " << setName
            << endl << endl;
    }

    // Mark all meshpoints on patch

    boolList vertOnPatch(mesh.nPoints(), false);

    forAll(meshPoints, pointI)
    {
        label meshPointI = meshPoints[pointI];

        vertOnPatch[meshPointI] = true;
    }


    // Mark cut edges.

    DynamicList<label> allCutEdges(pp.nEdges());

    DynamicList<scalar> allCutEdgeWeights(pp.nEdges());

    forAll(meshPoints, pointI)
    {
        label meshPointI = meshPoints[pointI];

        const labelList& pEdges = mesh.pointEdges()[meshPointI];

        forAll(pEdges, pEdgeI)
        {
            label edgeI = pEdges[pEdgeI];

            const edge& e = mesh.edges()[edgeI];

            label otherPointI = e.otherVertex(meshPointI);

            if (!vertOnPatch[otherPointI])
            {
                allCutEdges.append(edgeI);

                if (e.start() == meshPointI)
                {
                    allCutEdgeWeights.append(weight);
                }
                else
                {
                    allCutEdgeWeights.append(1 - weight);
                }
            }
        }
    }

    allCutEdges.shrink();
    allCutEdgeWeights.shrink();

    Info<< "Cutting:" << endl
        << "    cells:" << cutCells.size() << endl
        << "    edges:" << allCutEdges.size() << endl
        << endl;

    // Transfer DynamicLists to straight ones.
    scalarField cutEdgeWeights;
    cutEdgeWeights.transfer(allCutEdgeWeights);
    allCutEdgeWeights.clear();


    // Gets cuts across cells from cuts through edges.
    cellCuts cuts
    (
        mesh,
        cutCells.toc(),     // cells candidate for cutting
        labelList(0),       // cut vertices
        allCutEdges,        // cut edges
        cutEdgeWeights      // weight on cut edges
    );

    polyTopoChange meshMod(mesh);

    // Cutting engine
    meshCutter cutter(mesh);

    // Insert mesh refinement into polyTopoChange.
    cutter.setRefinement(cuts, meshMod);

    // Do all changes
    Info<< "Morphing ..." << endl;

    if (!overwrite)
    {
        runTime++;
    }

    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    // Update stored labels on meshCutter.
    cutter.updateMesh(morphMap());

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info << "Writing refined morphMesh to time " << runTime.timeName() << endl;

    mesh.write();

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
