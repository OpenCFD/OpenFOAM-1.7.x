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

Description
    'Stitches' a mesh.

    Takes a mesh and two patches and merges the faces on the two patches
    (if geometrically possible) so the faces become internal.

    Can do
    - 'perfect' match: faces and points on patches align exactly. Order might
    be different though.
    - 'integral' match: where the surfaces on both patches exactly
    match but the individual faces not
    - 'partial' match: where the non-overlapping part of the surface remains
    in the respective patch.

    Note : Is just a front-end to perfectInterface/slidingInterface.

    Comparable to running a meshModifier of the form
    (if masterPatch is called "M" and slavePatch "S"):

    couple
    {
        type                    slidingInterface;
        masterFaceZoneName      MSMasterZone
        slaveFaceZoneName       MSSlaveZone
        cutPointZoneName        MSCutPointZone
        cutFaceZoneName         MSCutFaceZone
        masterPatchName         M;
        slavePatchName          S;
        typeOfMatch             partial or integral
    }


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "ListOps.H"
#include "slidingInterface.H"
#include "perfectInterface.H"
#include "IOobjectList.H"
#include "ReadFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Checks whether patch present
void checkPatch(const polyBoundaryMesh& bMesh, const word& name)
{
    label patchI = bMesh.findPatchID(name);

    if (patchI == -1)
    {
        FatalErrorIn("checkPatch(const polyBoundaryMesh&, const word&)")
            << "Cannot find patch " << name << endl
            << "It should be present and of non-zero size" << endl
            << "Valid patches are " << bMesh.names()
            << exit(FatalError);
    }

    if (bMesh[patchI].empty())
    {
        FatalErrorIn("checkPatch(const polyBoundaryMesh&, const word&)")
            << "Patch " << name << " is present but zero size"
            << exit(FatalError);
    }
}


// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::noParallel();

    Foam::argList::validArgs.append("masterPatch");
    Foam::argList::validArgs.append("slavePatch");

    Foam::argList::validOptions.insert("partial", "");
    Foam::argList::validOptions.insert("perfect", "");

    Foam::argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();


    word masterPatchName(args.additionalArgs()[0]);
    word slavePatchName(args.additionalArgs()[1]);

    bool partialCover = args.optionFound("partial");
    bool perfectCover = args.optionFound("perfect");
    bool overwrite    = args.optionFound("overwrite");

    if (partialCover && perfectCover)
    {
        FatalErrorIn(args.executable())
            << "Cannot both supply partial and perfect." << endl
            << "Use perfect match option if the patches perfectly align"
            << " (both vertex positions and face centres)" << endl
            << exit(FatalError);
    }


    const word mergePatchName(masterPatchName + slavePatchName);
    const word cutZoneName(mergePatchName + "CutFaceZone");

    slidingInterface::typeOfMatch tom = slidingInterface::INTEGRAL;

    if (partialCover)
    {
        Info<< "Coupling partially overlapping patches "
            << masterPatchName << " and " << slavePatchName << nl
            << "Resulting internal faces will be in faceZone " << cutZoneName
            << nl
            << "Any uncovered faces will remain in their patch"
            << endl;

        tom = slidingInterface::PARTIAL;
    }
    else if (perfectCover)
    {
        Info<< "Coupling perfectly aligned patches "
            << masterPatchName << " and " << slavePatchName << nl
            << "Resulting (internal) faces will be in faceZone " << cutZoneName
            << nl << nl
            << "Note: both patches need to align perfectly." << nl
            << "Both the vertex"
            << " positions and the face centres need to align to within" << nl
            << "a tolerance given by the minimum edge length on the patch"
            << endl;
    }
    else
    {
        Info<< "Coupling patches " << masterPatchName << " and "
            << slavePatchName << nl
            << "Resulting (internal) faces will be in faceZone " << cutZoneName
            << nl << nl
            << "Note: the overall area covered by both patches should be"
            << " identical (\"integral\" interface)." << endl
            << "If this is not the case use the -partial option" << nl << endl;
    }

    // Check for non-empty master and slave patches
    checkPatch(mesh.boundaryMesh(), masterPatchName);
    checkPatch(mesh.boundaryMesh(), slavePatchName);

    // Create and add face zones and mesh modifiers

    // Master patch
    const polyPatch& masterPatch =
        mesh.boundaryMesh()
        [
            mesh.boundaryMesh().findPatchID(masterPatchName)
        ];

    // Make list of masterPatch faces
    labelList isf(masterPatch.size());

    forAll (isf, i)
    {
        isf[i] = masterPatch.start() + i;
    }

    polyTopoChanger stitcher(mesh);
    stitcher.setSize(1);

    DynamicList<pointZone*> pz;
    DynamicList<faceZone*> fz;
    DynamicList<cellZone*> cz;

    if (perfectCover)
    {
        // Add empty zone for resulting internal faces
        fz.append
        (
            new faceZone
            (
                cutZoneName,
                isf,
                boolList(masterPatch.size(), false),
                0,
                mesh.faceZones()
            )
        );

        // Note: make sure to add the zones BEFORE constructing polyMeshModifier
        // (since looks up various zones at construction time)
        Info << "Adding point and face zones" << endl;
        mesh.addZones(pz.shrink(), fz.shrink(), cz.shrink());

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
                masterPatchName,
                slavePatchName
            )
        );
    }
    else
    {
        pz.append
        (
            new pointZone
            (
                mergePatchName + "CutPointZone",
                labelList(0),
                0,
                mesh.pointZones()
            )
        );

        fz.append
        (
            new faceZone
            (
                mergePatchName + "MasterZone",
                isf,
                boolList(masterPatch.size(), false),
                0,
                mesh.faceZones()
            )
        );

        // Slave patch
        const polyPatch& slavePatch =
            mesh.boundaryMesh()
            [
                mesh.boundaryMesh().findPatchID(slavePatchName)
            ];

        labelList osf(slavePatch.size());

        forAll (osf, i)
        {
            osf[i] = slavePatch.start() + i;
        }

        fz.append
        (
            new faceZone
            (
                mergePatchName + "SlaveZone",
                osf,
                boolList(slavePatch.size(), false),
                1,
                mesh.faceZones()
            )
        );

        // Add empty zone for cut faces
        fz.append
        (
            new faceZone
            (
                cutZoneName,
                labelList(0),
                boolList(0, false),
                2,
                mesh.faceZones()
            )
        );


        // Note: make sure to add the zones BEFORE constructing polyMeshModifier
        // (since looks up various zones at construction time)
        Info << "Adding point and face zones" << endl;
        mesh.addZones(pz.shrink(), fz.shrink(), cz.shrink());

        // Add the sliding interface mesh modifier
        stitcher.set
        (
            0,
            new slidingInterface
            (
                "couple",
                0,
                stitcher,
                mergePatchName + "MasterZone",
                mergePatchName + "SlaveZone",
                mergePatchName + "CutPointZone",
                cutZoneName,
                masterPatchName,
                slavePatchName,
                tom,                    // integral or partial
                true                    // couple/decouple mode
            )
        );
    }


    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.timeName());

    // Read all current fvFields so they will get mapped
    Info<< "Reading all current volfields" << endl;
    PtrList<volScalarField> volScalarFields;
    ReadFields(mesh, objects, volScalarFields);

    PtrList<volVectorField> volVectorFields;
    ReadFields(mesh, objects, volVectorFields);

    PtrList<volSphericalTensorField> volSphericalTensorFields;
    ReadFields(mesh, objects, volSphericalTensorFields);

    PtrList<volSymmTensorField> volSymmTensorFields;
    ReadFields(mesh, objects, volSymmTensorFields);

    PtrList<volTensorField> volTensorFields;
    ReadFields(mesh, objects, volTensorFields);

    //- uncomment if you want to interpolate surface fields (usually bad idea)
    //Info<< "Reading all current surfaceFields" << endl;
    //PtrList<surfaceScalarField> surfaceScalarFields;
    //ReadFields(mesh, objects, surfaceScalarFields);
    //
    //PtrList<surfaceVectorField> surfaceVectorFields;
    //ReadFields(mesh, objects, surfaceVectorFields);
    //
    //PtrList<surfaceTensorField> surfaceTensorFields;
    //ReadFields(mesh, objects, surfaceTensorFields);

    if (!overwrite)
    {
        runTime++;
    }

    // Execute all polyMeshModifiers
    autoPtr<mapPolyMesh> morphMap = stitcher.changeMesh(true);

    mesh.movePoints(morphMap->preMotionPoints());

    // Write mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
        stitcher.instance() = oldInstance;
    }
    Info << nl << "Writing polyMesh to time " << runTime.timeName() << endl;

    IOstream::defaultPrecision(10);

    // Bypass runTime write (since only writes at outputTime)
    if
    (
       !runTime.objectRegistry::writeObject
        (
            runTime.writeFormat(),
            IOstream::currentVersion,
            runTime.writeCompression()
        )
    )
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    mesh.faceZones().write();
    mesh.pointZones().write();
    mesh.cellZones().write();

    // Write fields
    runTime.write();

    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //
