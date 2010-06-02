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

Application
    extrude2DMesh

Description
    Takes 2D mesh (all faces 2 points only, no front and back faces) and
    creates a 3D mesh by extruding with specified thickness.

Usage

    - extrude2DMesh thickness

    @param thickness \n
    Thickness (in metre) of slab.

Note
    Not sure about the walking of the faces to create the front and back faces.
    Tested on one .ccm file.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "extrude2DMesh.H"
#include "emptyPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("thickness");
    argList::validOptions.insert("overwrite", "");
#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    scalar thickness(readScalar(IStringStream(args.additionalArgs()[0])()));
    bool overwrite = args.optionFound("overwrite");


    // Check that mesh is 2D
    // ~~~~~~~~~~~~~~~~~~~~~

    const faceList& faces = mesh.faces();
    forAll(faces, faceI)
    {
        if (faces[faceI].size() != 2)
        {
            FatalErrorIn(args.executable())
                << "Face " << faceI << " size " << faces[faceI].size()
                << " is not of size 2 so mesh is not proper two-dimensional."
                << exit(FatalError);
        }
    }


    // Find extrude direction
    // ~~~~~~~~~~~~~~~~~~~~~~

    scalar minRange = GREAT;
    direction extrudeDir = 4;   //illegal value.

    for (direction dir = 0; dir < 3; dir++)
    {
        scalarField cmpts(mesh.points().component(dir));

        scalar range = max(cmpts)-min(cmpts);

        Info<< "Direction:" << dir << " range:" << range << endl;

        if (range < minRange)
        {
            minRange = range;
            extrudeDir = dir;
        }
    }

    Info<< "Extruding in direction " << extrudeDir
        << " with thickness " << thickness << nl
        << endl;



    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Add front and back patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    label frontPatchI = patches.findPatchID("frontAndBack");

    if (frontPatchI == -1)
    {
        // Add patch.
        List<polyPatch*> newPatches(patches.size()+1);

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            newPatches[patchI] = pp.clone
            (
                patches,
                newPatches.size(),
                pp.size(),
                pp.start()
            ).ptr();
        }

        frontPatchI = patches.size();

        newPatches[frontPatchI] = new emptyPolyPatch
        (
            "frontAndBack",
            0,
            mesh.nFaces(),
            frontPatchI,
            patches
        );

        Info<< "Adding empty patch " << newPatches[frontPatchI]->name()
            << " at index " << frontPatchI
            << " for front and back faces." << nl << endl;

        mesh.removeBoundary();
        mesh.addPatches(newPatches);
    }



    // Topo changes container. Initialise with number of patches.
    polyTopoChange meshMod(mesh.boundaryMesh().size());

    // Engine to extrude mesh
    extrude2DMesh extruder(mesh);

    // Insert changes into meshMod
    extruder.setRefinement
    (
        extrudeDir,
        thickness,
        frontPatchI,
        meshMod
    );

    // Create a mesh from topo changes.
    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    mesh.updateMesh(morphMap);

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    // Take over refinement levels and write to new time directory.
    Pout<< "Writing extruded mesh to time " << runTime.timeName() << nl
        << endl;

    mesh.write();

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
