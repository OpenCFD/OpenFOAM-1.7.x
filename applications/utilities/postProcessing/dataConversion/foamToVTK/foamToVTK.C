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

Application
    foamToVTK

Description
    Legacy VTK file format writer.

    - handles volScalar, volVector, pointScalar, pointVector, surfaceScalar
      fields.
    - mesh topo changes.
    - both ascii and binary.
    - single time step writing.
    - write subset only.
    - automatic decomposition of cells; polygons on boundary undecomposed since
      handled by vtk.

Usage

    - foamToVTK [OPTION]


    @param -ascii \n
    Write VTK data in ASCII format instead of binary.

    @param -mesh \<name\>\n
    Use a different mesh name (instead of -region)

    @param -fields \<fields\>\n
    Convert selected fields only. For example,
    @verbatim
         -fields "( p T U )"
    @endverbatim
    The quoting is required to avoid shell expansions and to pass the
    information as a single argument.

    @param -surfaceFields \n
    Write surfaceScalarFields (e.g., phi)

    @param -cellSet \<name\>\n
    @param -faceSet \<name\>\n
    @param -pointSet \<name\>\n
    Restrict conversion to the cellSet, faceSet or pointSet.

    @param -nearCellValue \n
    Output cell value on patches instead of patch value itself

    @param -noInternal \n
    Do not generate file for mesh, only for patches

    @param -noPointValues \n
    No pointFields

    @param -noFaceZones \n
    No faceZones

    @param -noLinks \n
    (in parallel) do not link processor files to master

    @param -allPatches \n
    Combine all patches into a single file

    @param -excludePatches \<patchNames\>\n
    Specify patches to exclude. For example,
    @verbatim
         -excludePatches "( inlet_1 inlet_2 )"
    @endverbatim
    The quoting is required to avoid shell expansions and to pass the
    information as a single argument.

    @param -useTimeName \n
    use the time index in the VTK file name instead of the time index

Note
    mesh subset is handled by vtkMesh. Slight inconsistency in
    interpolation: on the internal field it interpolates the whole volfield
    to the whole-mesh pointField and then selects only those values it
    needs for the subMesh (using the fvMeshSubset cellMap(), pointMap()
    functions). For the patches however it uses the
    fvMeshSubset.interpolate function to directly interpolate the
    whole-mesh values onto the subset patch.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "emptyPolyPatch.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"
#include "faceZoneMesh.H"
#include "Cloud.H"
#include "passiveParticle.H"

#include "vtkMesh.H"
#include "readFields.H"
#include "writeFuns.H"

#include "internalWriter.H"
#include "patchWriter.H"
#include "lagrangianWriter.H"

#include "writeFaceSet.H"
#include "writePointSet.H"
#include "writePatchGeom.H"
#include "writeSurfFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


static const label VTK_TETRA      = 10;
static const label VTK_PYRAMID    = 14;
static const label VTK_WEDGE      = 13;
static const label VTK_HEXAHEDRON = 12;


template<class GeoField>
void print(const char* msg, Ostream& os, const PtrList<GeoField>& flds)
{
    if (flds.size())
    {
        os << msg;
        forAll(flds, i)
        {
            os<< ' ' << flds[i].name();
        }
        os << endl;
    }
}


void print(Ostream& os, const wordList& flds)
{
    forAll(flds, i)
    {
        os<< ' ' << flds[i];
    }
    os << endl;
}


labelList getSelectedPatches
(
    const polyBoundaryMesh& patches,
    const HashSet<word>& excludePatches
)
{
    DynamicList<label> patchIDs(patches.size());

    Info<< "Combining patches:" << endl;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if
        (
            isA<emptyPolyPatch>(pp)
            || (Pstream::parRun() && isA<processorPolyPatch>(pp))
        )
        {
            Info<< "    discarding empty/processor patch " << patchI
                << " " << pp.name() << endl;
        }
        else if (!excludePatches.found(pp.name()))
        {
            patchIDs.append(patchI);
            Info<< "    patch " << patchI << " " << pp.name() << endl;
        }
    }
    return patchIDs.shrink();
}




// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "addRegionOption.H"

    argList::validOptions.insert("fields", "fields");
    argList::validOptions.insert("cellSet", "cellSet name");
    argList::validOptions.insert("faceSet", "faceSet name");
    argList::validOptions.insert("pointSet", "pointSet name");
    argList::validOptions.insert("ascii","");
    argList::validOptions.insert("surfaceFields","");
    argList::validOptions.insert("nearCellValue","");
    argList::validOptions.insert("noInternal","");
    argList::validOptions.insert("noPointValues","");
    argList::validOptions.insert("allPatches","");
    argList::validOptions.insert("excludePatches","patches to exclude");
    argList::validOptions.insert("noFaceZones","");
    argList::validOptions.insert("noLinks","");
    argList::validOptions.insert("useTimeName","");

#   include "setRootCase.H"
#   include "createTime.H"

    bool doWriteInternal = !args.optionFound("noInternal");
    bool doFaceZones     = !args.optionFound("noFaceZones");
    bool doLinks         = !args.optionFound("noLinks");
    bool binary          = !args.optionFound("ascii");
    bool useTimeName     = args.optionFound("useTimeName");

    if (binary && (sizeof(floatScalar) != 4 || sizeof(label) != 4))
    {
        FatalErrorIn(args.executable())
            << "floatScalar and/or label are not 4 bytes in size" << nl
            << "Hence cannot use binary VTK format. Please use -ascii"
            << exit(FatalError);
    }

    bool nearCellValue = args.optionFound("nearCellValue");

    if (nearCellValue)
    {
        WarningIn(args.executable())
            << "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }

    bool noPointValues = args.optionFound("noPointValues");

    if (noPointValues)
    {
        WarningIn(args.executable())
            << "Outputting cell values only" << nl << endl;
    }

    bool allPatches = args.optionFound("allPatches");

    HashSet<word> excludePatches;
    if (args.optionFound("excludePatches"))
    {
        args.optionLookup("excludePatches")() >> excludePatches;

        Info<< "Not including patches " << excludePatches << nl << endl;
    }

    word cellSetName;
    string vtkName;

    if (args.optionFound("cellSet"))
    {
        cellSetName = args.option("cellSet");
        vtkName = cellSetName;
    }
    else if (Pstream::parRun())
    {
        // Strip off leading casename, leaving just processor_DDD ending.
        vtkName = runTime.caseName();

        string::size_type i = vtkName.rfind("processor");

        if (i != string::npos)
        {
            vtkName = vtkName.substr(i);
        }
    }
    else
    {
        vtkName = runTime.caseName();
    }


    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"

    // VTK/ directory in the case
    fileName fvPath(runTime.path()/"VTK");
    // Directory of mesh (region0 gets filtered out)
    fileName regionPrefix = "";

    if (regionName != polyMesh::defaultRegion)
    {
        fvPath = fvPath/regionName;
        regionPrefix = regionName;
    }

    if (isDir(fvPath))
    {
        if
        (
            args.optionFound("time")
         || args.optionFound("latestTime")
         || cellSetName.size()
         || regionName != polyMesh::defaultRegion
        )
        {
            Info<< "Keeping old VTK files in " << fvPath << nl << endl;
        }
        else
        {
            Info<< "Deleting old VTK files in " << fvPath << nl << endl;

            rmDir(fvPath);
        }
    }

    mkDir(fvPath);


    // mesh wrapper; does subsetting and decomposition
    vtkMesh vMesh(mesh, cellSetName);


    // Scan for all possible lagrangian clouds
    HashSet<fileName> allCloudDirs;
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        fileNameList cloudDirs
        (
            readDir
            (
                runTime.timePath()/regionPrefix/cloud::prefix,
                fileName::DIRECTORY
            )
        );
        forAll(cloudDirs, i)
        {
            IOobjectList sprayObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudDirs[i]
            );

            IOobject* positionsPtr = sprayObjs.lookup("positions");

            if (positionsPtr)
            {
                if (allCloudDirs.insert(cloudDirs[i]))
                {
                    Info<< "At time: " << runTime.timeName()
                        << " detected cloud directory : " << cloudDirs[i]
                        << endl;
                }
            }
        }
    }


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time: " << runTime.timeName() << endl;

        word timeDesc = "";
        if (useTimeName)
        {
            timeDesc = runTime.timeName();
        }
        else
        {
            timeDesc = name(runTime.timeIndex());
        }

        // Check for new polyMesh/ and update mesh, fvMeshSubset and cell
        // decomposition.
        polyMesh::readUpdateState meshState = vMesh.readUpdate();

        const fvMesh& mesh = vMesh.mesh();

        if
        (
            meshState == polyMesh::TOPO_CHANGE
         || meshState == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            Info<< "    Read new mesh" << nl << endl;
        }


        // If faceSet: write faceSet only (as polydata)
        if (args.optionFound("faceSet"))
        {
            // Load the faceSet
            faceSet set(mesh, args.option("faceSet"));

            // Filename as if patch with same name.
            mkDir(fvPath/set.name());

            fileName patchFileName
            (
                fvPath/set.name()/set.name()
              + "_"
              + timeDesc
              + ".vtk"
            );

            Info<< "    FaceSet   : " << patchFileName << endl;

            writeFaceSet(binary, vMesh, set, patchFileName);

            continue;
        }
        // If pointSet: write pointSet only (as polydata)
        if (args.optionFound("pointSet"))
        {
            // Load the pointSet
            pointSet set(mesh, args.option("pointSet"));

            // Filename as if patch with same name.
            mkDir(fvPath/set.name());

            fileName patchFileName
            (
                fvPath/set.name()/set.name()
              + "_"
              + timeDesc
              + ".vtk"
            );

            Info<< "    pointSet   : " << patchFileName << endl;

            writePointSet(binary, vMesh, set, patchFileName);

            continue;
        }


        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        HashSet<word> selectedFields;
        if (args.optionFound("fields"))
        {
            args.optionLookup("fields")() >> selectedFields;
        }

        // Construct the vol fields (on the original mesh if subsetted)

        PtrList<volScalarField> vsf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vsf);
        print("    volScalarFields            :", Info, vsf);

        PtrList<volVectorField> vvf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vvf);
        print("    volVectorFields            :", Info, vvf);

        PtrList<volSphericalTensorField> vSpheretf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vSpheretf);
        print("    volSphericalTensorFields   :", Info, vSpheretf);

        PtrList<volSymmTensorField> vSymmtf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vSymmtf);
        print("    volSymmTensorFields        :", Info, vSymmtf);

        PtrList<volTensorField> vtf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vtf);
        print("    volTensorFields            :", Info, vtf);

        label nVolFields =
                vsf.size()
              + vvf.size()
              + vSpheretf.size()
              + vSymmtf.size()
              + vtf.size();


        // Construct pointMesh only if nessecary since constructs edge
        // addressing (expensive on polyhedral meshes)
        if (noPointValues)
        {
            Info<< "    pointScalarFields : switched off"
                << " (\"-noPointValues\" option)\n";
            Info<< "    pointVectorFields : switched off"
                << " (\"-noPointValues\" option)\n";
        }

        PtrList<pointScalarField> psf;
        PtrList<pointVectorField> pvf;
        PtrList<pointSphericalTensorField> pSpheretf;
        PtrList<pointSymmTensorField> pSymmtf;
        PtrList<pointTensorField> ptf;

        if (!noPointValues)
        {
            readFields
            (
                vMesh,
                pointMesh::New(vMesh.baseMesh()),
                objects,
                selectedFields,
                psf
            );
            print("    pointScalarFields          :", Info, psf);

            readFields
            (
                vMesh,
                pointMesh::New(vMesh.baseMesh()),
                objects,
                selectedFields,
                pvf
            );
            print("    pointVectorFields          :", Info, pvf);

            readFields
            (
                vMesh,
                pointMesh::New(vMesh.baseMesh()),
                objects,
                selectedFields,
                pSpheretf
            );
            print("    pointSphericalTensorFields :", Info, pSpheretf);

            readFields
            (
                vMesh,
                pointMesh::New(vMesh.baseMesh()),
                objects,
                selectedFields,
                pSymmtf
            );
            print("    pointSymmTensorFields      :", Info, pSymmtf);

            readFields
            (
                vMesh,
                pointMesh::New(vMesh.baseMesh()),
                objects,
                selectedFields,
                ptf
            );
            print("    pointTensorFields          :", Info, ptf);
        }
        Info<< endl;

        label nPointFields =
            psf.size()
          + pvf.size()
          + pSpheretf.size()
          + pSymmtf.size()
          + ptf.size();

        if (doWriteInternal)
        {
            //
            // Create file and write header
            //
            fileName vtkFileName
            (
                fvPath/vtkName
              + "_"
              + timeDesc
              + ".vtk"
            );

            Info<< "    Internal  : " << vtkFileName << endl;

            // Write mesh
            internalWriter writer(vMesh, binary, vtkFileName);

            // VolFields + cellID
            writeFuns::writeCellDataHeader
            (
                writer.os(),
                vMesh.nFieldCells(),
                1+nVolFields
            );

            // Write cellID field
            writer.writeCellIDs();

            // Write volFields
            writer.write(vsf);
            writer.write(vvf);
            writer.write(vSpheretf);
            writer.write(vSymmtf);
            writer.write(vtf);

            if (!noPointValues)
            {
                writeFuns::writePointDataHeader
                (
                    writer.os(),
                    vMesh.nFieldPoints(),
                    nVolFields+nPointFields
                );

                // pointFields
                writer.write(psf);
                writer.write(pvf);
                writer.write(pSpheretf);
                writer.write(pSymmtf);
                writer.write(ptf);

                // Interpolated volFields
                volPointInterpolation pInterp(mesh);
                writer.write(pInterp, vsf);
                writer.write(pInterp, vvf);
                writer.write(pInterp, vSpheretf);
                writer.write(pInterp, vSymmtf);
                writer.write(pInterp, vtf);
            }
        }

        //---------------------------------------------------------------------
        //
        // Write surface fields
        //
        //---------------------------------------------------------------------

        if (args.optionFound("surfaceFields"))
        {
            PtrList<surfaceScalarField> ssf;
            readFields
            (
                vMesh,
                vMesh.baseMesh(),
                objects,
                selectedFields,
                ssf
            );
            print("    surfScalarFields  :", Info, ssf);

            PtrList<surfaceVectorField> svf;
            readFields
            (
                vMesh,
                vMesh.baseMesh(),
                objects,
                selectedFields,
                svf
            );
            print("    surfVectorFields  :", Info, svf);

            if (ssf.size() + svf.size() > 0)
            {
                // Rework the scalar fields into vectorfields.
                label sz = svf.size();

                svf.setSize(sz+ssf.size());

                surfaceVectorField n(mesh.Sf()/mesh.magSf());

                forAll(ssf, i)
                {
                    svf.set(sz+i, ssf[i]*n);
                    svf[sz+i].rename(ssf[i].name());
                }
                ssf.clear();

                mkDir(fvPath / "surfaceFields");

                fileName surfFileName
                (
                    fvPath
                   /"surfaceFields"
                   /"surfaceFields"
                   + "_"
                   + timeDesc
                   + ".vtk"
                );

                writeSurfFields
                (
                    binary,
                    vMesh,
                    surfFileName,
                    svf
                );
            }
        }


        //---------------------------------------------------------------------
        //
        // Write patches (POLYDATA file, one for each patch)
        //
        //---------------------------------------------------------------------

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        if (allPatches)
        {
            mkDir(fvPath/"allPatches");

            fileName patchFileName;

            if (vMesh.useSubMesh())
            {
                patchFileName =
                    fvPath/"allPatches"/cellSetName
                  + "_"
                  + timeDesc
                  + ".vtk";
            }
            else
            {
                patchFileName =
                    fvPath/"allPatches"/"allPatches"
                  + "_"
                  + timeDesc
                  + ".vtk";
            }

            Info<< "    Combined patches     : " << patchFileName << endl;

            patchWriter writer
            (
                vMesh,
                binary,
                nearCellValue,
                patchFileName,
                getSelectedPatches(patches, excludePatches)
            );

            // VolFields + patchID
            writeFuns::writeCellDataHeader
            (
                writer.os(),
                writer.nFaces(),
                1+nVolFields
            );

            // Write patchID field
            writer.writePatchIDs();

            // Write volFields
            writer.write(vsf);
            writer.write(vvf);
            writer.write(vSpheretf);
            writer.write(vSymmtf);
            writer.write(vtf);

            if (!noPointValues)
            {
                writeFuns::writePointDataHeader
                (
                    writer.os(),
                    writer.nPoints(),
                    nPointFields
                );

                // Write pointFields
                writer.write(psf);
                writer.write(pvf);
                writer.write(pSpheretf);
                writer.write(pSymmtf);
                writer.write(ptf);

                // no interpolated volFields since I cannot be bothered to
                // create the patchInterpolation for all subpatches.
            }
        }
        else
        {
            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (!excludePatches.found(pp.name()))
                {
                    mkDir(fvPath/pp.name());

                    fileName patchFileName;

                    if (vMesh.useSubMesh())
                    {
                        patchFileName =
                            fvPath/pp.name()/cellSetName
                          + "_"
                          + timeDesc
                          + ".vtk";
                    }
                    else
                    {
                        patchFileName =
                            fvPath/pp.name()/pp.name()
                          + "_"
                          + timeDesc
                          + ".vtk";
                    }

                    Info<< "    Patch     : " << patchFileName << endl;

                    patchWriter writer
                    (
                        vMesh,
                        binary,
                        nearCellValue,
                        patchFileName,
                        labelList(1, patchI)
                    );

                    if (!isA<emptyPolyPatch>(pp))
                    {
                        // VolFields + patchID
                        writeFuns::writeCellDataHeader
                        (
                            writer.os(),
                            writer.nFaces(),
                            1+nVolFields
                        );

                        // Write patchID field
                        writer.writePatchIDs();

                        // Write volFields
                        writer.write(vsf);
                        writer.write(vvf);
                        writer.write(vSpheretf);
                        writer.write(vSymmtf);
                        writer.write(vtf);

                        if (!noPointValues)
                        {
                            writeFuns::writePointDataHeader
                            (
                                writer.os(),
                                writer.nPoints(),
                                nVolFields
                              + nPointFields
                            );

                            // Write pointFields
                            writer.write(psf);
                            writer.write(pvf);
                            writer.write(pSpheretf);
                            writer.write(pSymmtf);
                            writer.write(ptf);

                            PrimitivePatchInterpolation<primitivePatch> pInter
                            (
                                pp
                            );

                            // Write interpolated volFields
                            writer.write(pInter, vsf);
                            writer.write(pInter, vvf);
                            writer.write(pInter, vSpheretf);
                            writer.write(pInter, vSymmtf);
                            writer.write(pInter, vtf);
                        }
                    }
                }
            }
        }

        //---------------------------------------------------------------------
        //
        // Write faceZones (POLYDATA file, one for each zone)
        //
        //---------------------------------------------------------------------

        if (doFaceZones)
        {
            const faceZoneMesh& zones = mesh.faceZones();

            forAll(zones, zoneI)
            {
                const faceZone& pp = zones[zoneI];

                mkDir(fvPath/pp.name());

                fileName patchFileName;

                if (vMesh.useSubMesh())
                {
                    patchFileName =
                        fvPath/pp.name()/cellSetName
                      + "_"
                      + timeDesc
                      + ".vtk";
                }
                else
                {
                    patchFileName =
                        fvPath/pp.name()/pp.name()
                      + "_"
                      + timeDesc
                      + ".vtk";
                }

                Info<< "    FaceZone  : " << patchFileName << endl;

                std::ofstream str(patchFileName.c_str());

                writeFuns::writeHeader(str, binary, pp.name());
                str << "DATASET POLYDATA" << std::endl;

                writePatchGeom
                (
                    binary,
                    pp().localFaces(),
                    pp().localPoints(),
                    str
                );
            }
        }



        //---------------------------------------------------------------------
        //
        // Write lagrangian data
        //
        //---------------------------------------------------------------------

        forAllConstIter(HashSet<fileName>, allCloudDirs, iter)
        {
            const fileName& cloudName = iter.key();

            // Always create the cloud directory.
            mkDir(fvPath/cloud::prefix/cloudName);

            fileName lagrFileName
            (
                fvPath/cloud::prefix/cloudName/cloudName
              + "_" + timeDesc + ".vtk"
            );

            Info<< "    Lagrangian: " << lagrFileName << endl;


            IOobjectList sprayObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudName
            );

            IOobject* positionsPtr = sprayObjs.lookup("positions");

            if (positionsPtr)
            {
                wordList labelNames(sprayObjs.names(labelIOField::typeName));
                Info<< "        labels            :";
                print(Info, labelNames);

                wordList scalarNames(sprayObjs.names(scalarIOField::typeName));
                Info<< "        scalars           :";
                print(Info, scalarNames);

                wordList vectorNames(sprayObjs.names(vectorIOField::typeName));
                Info<< "        vectors           :";
                print(Info, vectorNames);

                wordList sphereNames
                (
                    sprayObjs.names
                    (
                        sphericalTensorIOField::typeName
                    )
                );
                Info<< "        spherical tensors :";
                print(Info, sphereNames);

                wordList symmNames
                (
                    sprayObjs.names
                    (
                        symmTensorIOField::typeName
                    )
                );
                Info<< "        symm tensors      :";
                print(Info, symmNames);

                wordList tensorNames(sprayObjs.names(tensorIOField::typeName));
                Info<< "        tensors           :";
                print(Info, tensorNames);

                lagrangianWriter writer
                (
                    vMesh,
                    binary,
                    lagrFileName,
                    cloudName,
                    false
                );

                // Write number of fields
                writer.writeParcelHeader
                (
                    labelNames.size()
                  + scalarNames.size()
                  + vectorNames.size()
                  + sphereNames.size()
                  + symmNames.size()
                  + tensorNames.size()
                );

                // Fields
                writer.writeIOField<label>(labelNames);
                writer.writeIOField<scalar>(scalarNames);
                writer.writeIOField<vector>(vectorNames);
                writer.writeIOField<sphericalTensor>(sphereNames);
                writer.writeIOField<symmTensor>(symmNames);
                writer.writeIOField<tensor>(tensorNames);
            }
            else
            {
                lagrangianWriter writer
                (
                    vMesh,
                    binary,
                    lagrFileName,
                    cloudName,
                    true
                );

                // Write number of fields
                writer.writeParcelHeader(0);
            }
        }
    }


    //---------------------------------------------------------------------
    //
    // Link parallel outputs back to undecomposed case for ease of loading
    //
    //---------------------------------------------------------------------

    if (Pstream::parRun() && doLinks)
    {
        mkDir(runTime.path()/".."/"VTK");
        chDir(runTime.path()/".."/"VTK");

        Info<< "Linking all processor files to " << runTime.path()/".."/"VTK"
            << endl;

        // Get list of vtk files
        fileName procVTK
        (
            fileName("..")
           /"processor" + name(Pstream::myProcNo())
           /"VTK"
        );

        fileNameList dirs(readDir(procVTK, fileName::DIRECTORY));
        label sz = dirs.size();
        dirs.setSize(sz+1);
        dirs[sz] = ".";

        forAll(dirs, i)
        {
            fileNameList subFiles(readDir(procVTK/dirs[i], fileName::FILE));

            forAll(subFiles, j)
            {
                fileName procFile(procVTK/dirs[i]/subFiles[j]);

                if (exists(procFile))
                {
                    string cmd
                    (
                        "ln -s "
                      + procFile
                      + " "
                      + "processor"
                      + name(Pstream::myProcNo())
                      + "_"
                      + procFile.name()
                    );
                    system(cmd.c_str());
                }
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
