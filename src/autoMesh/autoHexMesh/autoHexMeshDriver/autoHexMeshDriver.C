/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*----------------------------------------------------------------------------*/

#include "autoHexMeshDriver.H"
#include "fvMesh.H"
#include "Time.H"
#include "boundBox.H"
#include "wallPolyPatch.H"
#include "cellSet.H"
#include "syncTools.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "layerParameters.H"
#include "autoRefineDriver.H"
#include "autoSnapDriver.H"
#include "autoLayerDriver.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(autoHexMeshDriver, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Check writing tolerance before doing any serious work
Foam::scalar Foam::autoHexMeshDriver::getMergeDistance(const scalar mergeTol)
 const
{
    const boundBox& meshBb = mesh_.bounds();
    scalar mergeDist = mergeTol * meshBb.mag();
    scalar writeTol = std::pow
    (
        scalar(10.0),
       -scalar(IOstream::defaultPrecision())
    );

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    if (mesh_.time().writeFormat() == IOstream::ASCII && mergeTol < writeTol)
    {
        FatalErrorIn("autoHexMeshDriver::getMergeDistance(const scalar) const")
            << "Your current settings specify ASCII writing with "
            << IOstream::defaultPrecision() << " digits precision." << endl
            << "Your merging tolerance (" << mergeTol << ") is finer than this."
            << endl
            << "Please change your writeFormat to binary"
            << " or increase the writePrecision" << endl
            << "or adjust the merge tolerance (-mergeTol)."
            << exit(FatalError);
    }

    return mergeDist;
}


//// Specifically orient using a calculated point outside
//void Foam::autoHexMeshDriver::orientOutside
//(
//    PtrList<searchableSurface>& shells
//)
//{
//    // Determine outside point.
//    boundBox overallBb = boundBox::invertedBox;
//
//    bool hasSurface = false;
//
//    forAll(shells, shellI)
//    {
//        if (isA<triSurfaceMesh>(shells[shellI]))
//        {
//            const triSurfaceMesh& shell =
//                refCast<const triSurfaceMesh>(shells[shellI]);
//
//            hasSurface = true;
//
//            boundBox shellBb(shell.localPoints(), false);
//
//            overallBb.min() = min(overallBb.min(), shellBb.min());
//            overallBb.max() = max(overallBb.max(), shellBb.max());
//        }
//    }
//
//    if (hasSurface)
//    {
//        const point outsidePt = 2 * overallBb.span();
//
//        //Info<< "Using point " << outsidePt << " to orient shells" << endl;
//
//        forAll(shells, shellI)
//        {
//            if (isA<triSurfaceMesh>(shells[shellI]))
//            {
//                triSurfaceMesh& shell =
//                  refCast<triSurfaceMesh>(shells[shellI]);
//
//                if (!refinementSurfaces::isSurfaceClosed(shell))
//                {
//                    FatalErrorIn("orientOutside(PtrList<searchableSurface>&)")
//                        << "Refinement shell "
//                        << shell.searchableSurface::name()
//                        << " is not closed." << exit(FatalError);
//                }
//
//                refinementSurfaces::orientSurface(outsidePt, shell);
//            }
//        }
//    }
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::autoHexMeshDriver::autoHexMeshDriver
(
    fvMesh& mesh,
    const bool overwrite,
    const dictionary& dict,
    const dictionary& decomposeDict
)
:
    mesh_(mesh),
    dict_(dict),
    debug_(readLabel(dict_.lookup("debug"))),
    mergeDist_(getMergeDistance(readScalar(dict_.lookup("mergeTolerance"))))
{
    if (debug_ > 0)
    {
        meshRefinement::debug = debug_;
        autoHexMeshDriver::debug = debug_;
        autoRefineDriver::debug = debug;
        autoSnapDriver::debug = debug;
        autoLayerDriver::debug = debug;
    }

    refinementParameters refineParams(dict, 1);

    Info<< "Overall cell limit                         : "
        << refineParams.maxGlobalCells() << endl;
    Info<< "Per processor cell limit                   : "
        << refineParams.maxLocalCells() << endl;
    Info<< "Minimum number of cells to refine          : "
        << refineParams.minRefineCells() << endl;
    Info<< "Curvature                                  : "
        << refineParams.curvature() << nl << endl;
    Info<< "Layers between different refinement levels : "
        << refineParams.nBufferLayers() << endl;

    PtrList<dictionary> shellDicts(dict_.lookup("refinementShells"));

    PtrList<dictionary> surfaceDicts(dict_.lookup("surfaces"));


    // Read geometry
    // ~~~~~~~~~~~~~

    {
        Info<< "Reading all geometry." << endl;

        // Construct dictionary with all shells and all refinement surfaces
        dictionary geometryDict;

        forAll(shellDicts, shellI)
        {
            dictionary shellDict = shellDicts[shellI];
            const word name(shellDict.lookup("name"));
            shellDict.remove("name");
            shellDict.remove("level");
            shellDict.remove("refineInside");
            geometryDict.add(name, shellDict);
        }

        forAll(surfaceDicts, surfI)
        {
            dictionary surfDict = surfaceDicts[surfI];
            const word name(string::validate<word>(surfDict.lookup("file")));
            surfDict.remove("file");
            surfDict.remove("regions");
            if (!surfDict.found("name"))
            {
                surfDict.add("name", name);
            }
            surfDict.add("type", triSurfaceMesh::typeName);
            geometryDict.add(name, surfDict);
        }

        allGeometryPtr_.reset
        (
            new searchableSurfaces
            (
                IOobject
                (
                    "abc",                                      // dummy name
                    //mesh_.time().findInstance("triSurface", word::null),
                                                                // instance
                    mesh_.time().constant(),                    // instance
                    "triSurface",                               // local
                    mesh_.time(),                               // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                geometryDict
            )
        );

        Info<< "Read geometry in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;
    }


    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    {
        Info<< "Reading surfaces and constructing search trees." << endl;

        surfacesPtr_.reset
        (
            new refinementSurfaces
            (
                allGeometryPtr_(),
                surfaceDicts
            )
        );
        Info<< "Read surfaces in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;
    }

    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        Info<< "Reading refinement shells." << endl;
        shellsPtr_.reset
        (
            new shellSurfaces
            (
                allGeometryPtr_(),
                shellDicts
            )
        );
        Info<< "Read refinement shells in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        //// Orient shell surfaces before any searching is done.
        //Info<< "Orienting triSurface shells so point far away is outside."
        //    << endl;
        //orientOutside(shells_);
        //Info<< "Oriented shells in = "
        //    << mesh_.time().cpuTimeIncrement() << " s" << endl;

        Info<< "Setting refinement level of surface to be consistent"
            << " with shells." << endl;
        surfacesPtr_().setMinLevelFields(shells());
        Info<< "Checked shell refinement in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;
    }


    // Check faceZones are synchronised
    meshRefinement::checkCoupledFaceZones(mesh_);


    // Refinement engine
    // ~~~~~~~~~~~~~~~~~

    {
        Info<< nl
            << "Determining initial surface intersections" << nl
            << "-----------------------------------------" << nl
            << endl;

        // Main refinement engine
        meshRefinerPtr_.reset
        (
            new meshRefinement
            (
                mesh,
                mergeDist_,         // tolerance used in sorting coordinates
                overwrite,
                surfaces(),
                shells()
            )
        );
        Info<< "Calculated surface intersections in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // Some stats
        meshRefinerPtr_().printMeshInfo(debug_, "Initial mesh");

        meshRefinerPtr_().write
        (
            debug_&meshRefinement::OBJINTERSECTIONS,
            mesh_.time().path()/meshRefinerPtr_().timeName()
        );
    }


    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        Info<< nl
            << "Adding patches for surface regions" << nl
            << "----------------------------------" << nl
            << endl;

        // From global region number to mesh patch.
        globalToPatch_.setSize(surfaces().nRegions(), -1);

        Info<< "Patch\tRegion" << nl
            << "-----\t------"
            << endl;

        const labelList& surfaceGeometry = surfaces().surfaces();
        forAll(surfaceGeometry, surfI)
        {
            label geomI = surfaceGeometry[surfI];

            const wordList& regNames = allGeometryPtr_().regionNames()[geomI];

            Info<< surfaces().names()[surfI] << ':' << nl << nl;

            forAll(regNames, i)
            {
                label patchI = meshRefinerPtr_().addMeshedPatch
                (
                    regNames[i],
                    wallPolyPatch::typeName
                );

                Info<< patchI << '\t' << regNames[i] << nl;

                globalToPatch_[surfaces().globalRegion(surfI, i)] = patchI;
            }

            Info<< nl;
        }
        Info<< "Added patches in = "
            << mesh_.time().cpuTimeIncrement() << " s" << nl << endl;
    }


    //// Add cyclics for any named faceZones
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //// (these cyclics are used later on to temporarily put the faceZones
    ////  in when snapping)
    //
    //labelList namedSurfaces(surfaces().getNamedSurfaces());
    //if (namedSurfaces.size())
    //{
    //    Info<< nl
    //        << "Introducing cyclics for faceZones" << nl
    //        << "---------------------------------" << nl
    //        << endl;
    //
    //    // From surface to cyclic patch
    //    surfaceToCyclicPatch_.setSize(surfaces().size(), -1);
    //
    //    Info<< "Patch\tZone" << nl
    //        << "----\t-----"
    //        << endl;
    //
    //    forAll(namedSurfaces, i)
    //    {
    //        label surfI = namedSurfaces[i];
    //
    //        surfaceToCyclicPatch_[surfI] = meshRefinement::addPatch
    //        (
    //            mesh,
    //            surfaces().faceZoneNames()[surfI],
    //            cyclicPolyPatch::typeName
    //        );
    //
    //        Info<< surfaceToCyclicPatch_[surfI] << '\t'
    //            << surfaces().faceZoneNames()[surfI] << nl << endl;
    //    }
    //    Info<< "Added cyclic patches in = "
    //        << mesh_.time().cpuTimeIncrement() << " s" << endl;
    //}


    // Parallel
    // ~~~~~~~~

    {
        // Decomposition
        decomposerPtr_ = decompositionMethod::New
        (
            decomposeDict,
            mesh_
        );
        decompositionMethod& decomposer = decomposerPtr_();


        if (Pstream::parRun() && !decomposer.parallelAware())
        {
            FatalErrorIn("autoHexMeshDriver::autoHexMeshDriver"
                "(const IOobject&, fvMesh&)")
                << "You have selected decomposition method "
                << decomposer.typeName
                << " which is not parallel aware." << endl
                << "Please select one that is (parMetis, hierarchical)"
                << exit(FatalError);
        }

        // Mesh distribution engine (uses tolerance to reconstruct meshes)
        distributorPtr_.reset(new fvMeshDistribute(mesh_, mergeDist_));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::autoHexMeshDriver::writeMesh(const string& msg) const
{
    const meshRefinement& meshRefiner = meshRefinerPtr_();

    meshRefiner.printMeshInfo(debug_, msg);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    meshRefiner.write(meshRefinement::MESH|meshRefinement::SCALARLEVELS, "");
    if (debug_ & meshRefinement::OBJINTERSECTIONS)
    {
        meshRefiner.write
        (
            meshRefinement::OBJINTERSECTIONS,
            mesh_.time().path()/meshRefiner.timeName()
        );
    }
    Info<< "Written mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s." << endl;
}


void Foam::autoHexMeshDriver::doMesh()
{
    Switch wantRefine(dict_.lookup("doRefine"));
    Switch wantSnap(dict_.lookup("doSnap"));
    Switch wantLayers(dict_.lookup("doLayers"));

    Info<< "Do refinement : " << wantRefine << nl
        << "Do snapping   : " << wantSnap << nl
        << "Do layers     : " << wantLayers << nl
        << endl;

    if (wantRefine)
    {
        const dictionary& motionDict = dict_.subDict("motionDict");

        autoRefineDriver refineDriver
        (
            meshRefinerPtr_(),
            decomposerPtr_(),
            distributorPtr_(),
            globalToPatch_
        );

        // Get all the refinement specific params
        refinementParameters refineParams(dict_, 1);

        refineDriver.doRefine(dict_, refineParams, wantSnap, motionDict);

        // Write mesh
        writeMesh("Refined mesh");
    }

    if (wantSnap)
    {
        const dictionary& snapDict = dict_.subDict("snapDict");
        const dictionary& motionDict = dict_.subDict("motionDict");

        autoSnapDriver snapDriver
        (
            meshRefinerPtr_(),
            globalToPatch_
        );

        // Get all the snapping specific params
        snapParameters snapParams(snapDict, 1);

        snapDriver.doSnap(snapDict, motionDict, snapParams);

        // Write mesh.
        writeMesh("Snapped mesh");
    }

    if (wantLayers)
    {
        const dictionary& motionDict = dict_.subDict("motionDict");
        const dictionary& shrinkDict = dict_.subDict("shrinkDict");
        PtrList<dictionary> surfaceDicts(dict_.lookup("surfaces"));

        autoLayerDriver layerDriver(meshRefinerPtr_());

        // Get all the layer specific params
        layerParameters layerParams
        (
            surfaceDicts,
            surfacesPtr_(),
            globalToPatch_,
            shrinkDict,
            mesh_.boundaryMesh()
        );

        layerDriver.doLayers
        (
            shrinkDict,
            motionDict,
            layerParams,
            true,                   // pre-balance
            decomposerPtr_(),
            distributorPtr_()
        );

        // Write mesh.
        writeMesh("Layer mesh");
    }
}


// ************************************************************************* //
