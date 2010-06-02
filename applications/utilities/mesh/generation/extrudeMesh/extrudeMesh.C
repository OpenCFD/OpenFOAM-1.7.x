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
    Extrude mesh from existing patch (by default outwards facing normals;
    optional flips faces) or from patch read from file.

    Note: Merges close points so be careful.

    Type of extrusion prescribed by run-time selectable model.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dimensionedTypes.H"
#include "IFstream.H"
#include "faceMesh.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "edgeCollapser.H"
#include "mathematicalConstants.H"
#include "globalMeshData.H"
#include "perfectInterface.H"

#include "extrudedMesh.H"
#include "extrudeModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTimeExtruded.H"

    autoPtr<extrudedMesh> meshPtr(NULL);

    IOdictionary dict
    (
        IOobject
        (
            "extrudeProperties",
            runTimeExtruded.constant(),
            runTimeExtruded,
            IOobject::MUST_READ
        )
    );

    autoPtr<extrudeModel> model(extrudeModel::New(dict));

    const word sourceType(dict.lookup("constructFrom"));

    autoPtr<faceMesh> fMesh;

    if (sourceType == "patch")
    {
        fileName sourceCasePath(dict.lookup("sourceCase"));
        sourceCasePath.expand();
        fileName sourceRootDir = sourceCasePath.path();
        fileName sourceCaseDir = sourceCasePath.name();
        word patchName(dict.lookup("sourcePatch"));

        Info<< "Extruding patch " << patchName
            << " on mesh " << sourceCasePath << nl
            << endl;

        Time runTime
        (
            Time::controlDictName,
            sourceRootDir,
            sourceCaseDir
        );
        #include "createPolyMesh.H"

        label patchID = mesh.boundaryMesh().findPatchID(patchName);

        if (patchID == -1)
        {
            FatalErrorIn(args.executable())
                << "Cannot find patch " << patchName
                << " in the source mesh.\n"
                << "Valid patch names are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchID];
        fMesh.reset(new faceMesh(pp.localFaces(), pp.localPoints()));

        {
            fileName surfName(runTime.path()/patchName + ".sMesh");
            Info<< "Writing patch as surfaceMesh to "
                << surfName << nl << endl;
            OFstream os(surfName);
            os << fMesh() << nl;
        }
    }
    else if (sourceType == "surface")
    {
        // Read from surface
        fileName surfName(dict.lookup("surface"));

        Info<< "Extruding surfaceMesh read from file " << surfName << nl
            << endl;

        IFstream is(surfName);

        fMesh.reset(new faceMesh(is));

        Info<< "Read patch from file " << surfName << nl
            << endl;
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Illegal 'constructFrom' specification. Should either be "
            << "patch or surface." << exit(FatalError);
    }

    Switch flipNormals(dict.lookup("flipNormals"));

    if (flipNormals)
    {
        Info<< "Flipping faces." << nl << endl;

        faceList faces(fMesh().size());
        forAll(faces, i)
        {
            faces[i] = fMesh()[i].reverseFace();
        }
        fMesh.reset(new faceMesh(faces, fMesh().localPoints()));
    }


    Info<< "Extruding patch with :" << nl
            << "    points     : " << fMesh().points().size() << nl
            << "    faces      : " << fMesh().size() << nl
            << "    normals[0] : " << fMesh().faceNormals()[0]
            << nl
            << endl;

    extrudedMesh mesh
    (
        IOobject
        (
            extrudedMesh::defaultRegion,
            runTimeExtruded.constant(),
            runTimeExtruded
        ),
        fMesh(),
        model()
    );


    const boundBox& bb = mesh.globalData().bb();
    const vector span = bb.span();
    const scalar mergeDim = 1E-4 * bb.minDim();

    Info<< "Mesh bounding box : " << bb << nl
        << "        with span : " << span << nl
        << "Merge distance    : " << mergeDim << nl
        << endl;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const label origPatchID = patches.findPatchID("originalPatch");
    const label otherPatchID = patches.findPatchID("otherSide");

    if (origPatchID == -1 || otherPatchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Cannot find patch originalPatch or otherSide." << nl
            << "Valid patches are " << patches.names() << exit(FatalError);
    }

    // Collapse edges
    // ~~~~~~~~~~~~~~

    {
        Info<< "Collapsing edges < " << mergeDim << " ..." << nl << endl;

        // Edge collapsing engine
        edgeCollapser collapser(mesh);

        const edgeList& edges = mesh.edges();
        const pointField& points = mesh.points();

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            scalar d = e.mag(points);

            if (d < mergeDim)
            {
                Info<< "Merging edge " << e << " since length " << d
                    << " << " << mergeDim << nl;

                // Collapse edge to e[0]
                collapser.collapseEdge(edgeI, e[0]);
            }
        }

        // Topo change container
        polyTopoChange meshMod(mesh);
        // Put all modifications into meshMod
        bool anyChange = collapser.setRefinement(meshMod);

        if (anyChange)
        {
            // Construct new mesh from polyTopoChange.
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

            // Update fields
            mesh.updateMesh(map);

            // Move mesh (if inflation used)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
        }
    }


    // Merging front and back patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Switch mergeFaces(dict.lookup("mergeFaces"));
    if (mergeFaces)
    {
        Info<< "Assuming full 360 degree axisymmetric case;"
            << " stitching faces on patches "
            << patches[origPatchID].name() << " and "
            << patches[otherPatchID].name() << " together ..." << nl << endl;

        polyTopoChanger stitcher(mesh);
        stitcher.setSize(1);

        // Make list of masterPatch faces
        labelList isf(patches[origPatchID].size());

        forAll (isf, i)
        {
            isf[i] = patches[origPatchID].start() + i;
        }

        const word cutZoneName("originalCutFaceZone");

        List<faceZone*> fz
        (
            1,
            new faceZone
            (
                cutZoneName,
                isf,
                boolList(isf.size(), false),
                0,
                mesh.faceZones()
            )
        );

        mesh.addZones(List<pointZone*>(0), fz, List<cellZone*>(0));

        // Add the perfect interface mesh modifier
        stitcher.set
        (
            0,
            new perfectInterface
            (
                "couple",
                0,
                stitcher,
                cutZoneName,
                patches[origPatchID].name(),
                patches[otherPatchID].name()
            )
        );

        // Execute all polyMeshModifiers
        autoPtr<mapPolyMesh> morphMap = stitcher.changeMesh(true);

        mesh.movePoints(morphMap->preMotionPoints());
    }

    if (!mesh.write())
    {
        FatalErrorIn(args.executable()) << "Failed writing mesh"
            << exit(FatalError);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
