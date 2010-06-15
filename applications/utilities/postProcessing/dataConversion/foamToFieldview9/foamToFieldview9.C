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

Description
    Write out the OpenFOAM mesh in Version 3.0 Fieldview-UNS format (binary).

    See Fieldview Release 9 Reference Manual - Appendix D
    (Unstructured Data Format)
    Borrows various from uns/write_binary_uns.c from FieldView dist.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "scalarIOField.H"
#include "volPointInterpolation.H"
#include "wallFvPatch.H"
#include "symmetryFvPatch.H"

#include "Cloud.H"
#include "passiveParticle.H"

#include "IOobjectList.H"
#include "boolList.H"
#include "stringList.H"
#include "cellModeller.H"

#include "floatScalar.H"
#include "calcFaceAddressing.H"
#include "writeFunctions.H"
#include "fieldviewTopology.H"

#include <fstream>

#include "fv_reader_tags.h"

extern "C"
{
    unsigned int fv_encode_elem_header(int elem_type, int wall_info[]);
}

using namespace Foam;

typedef Field<floatScalar> floatField;


static HashTable<word> FieldviewNames;


static word getFieldViewName(const word& foamName)
{
    if (FieldviewNames.found(foamName))
    {
        return FieldviewNames[foamName];
    }
    else
    {
        return foamName;
    }
}


static void printNames(const wordList& names, Ostream& os)
{
    forAll(names, fieldI)
    {
        Info<< " " << names[fieldI] << '/' << getFieldViewName(names[fieldI]);
    }
}


// Count number of vertices used by celli
static label countVerts(const primitiveMesh& mesh, const label celli)
{
    const cell& cll = mesh.cells()[celli];

    // Count number of vertices used
    labelHashSet usedVerts(10*cll.size());

    forAll(cll, cellFacei)
    {
        const face& f = mesh.faces()[cll[cellFacei]];

        forAll(f, fp)
        {
            if (!usedVerts.found(f[fp]))
            {
                usedVerts.insert(f[fp]);
            }
        }
    }
    return usedVerts.toc().size();
}


static void writeFaceData
(
    const polyMesh& mesh,
    const fieldviewTopology& topo,
    const label patchI,
    const scalarField& patchField,
    const bool writePolyFaces,
    std::ofstream& fvFile
)
{
    const polyPatch& pp = mesh.boundaryMesh()[patchI];

    // Start off with dummy value.
    if (writePolyFaces)
    {
        floatField fField(topo.nPolyFaces()[patchI], 0.0);

        // Fill selected faces with field values
        label polyFaceI = 0;
        forAll(patchField, faceI)
        {
            if (pp[faceI].size() > 4)
            {
                fField[polyFaceI++] = float(patchField[faceI]);
            }
        }

        fvFile.write
        (
            reinterpret_cast<char*>(fField.begin()), fField.size()*sizeof(float)
        );
    }
    else
    {
        floatField fField(pp.size() - topo.nPolyFaces()[patchI], 0.0);

        // Fill selected faces with field values
        label quadFaceI = 0;
        forAll(patchField, faceI)
        {
            if (pp[faceI].size() <= 4)
            {
                fField[quadFaceI++] = float(patchField[faceI]);
            }
        }

        fvFile.write
        (
            reinterpret_cast<char*>(fField.begin()), fField.size()*sizeof(float)
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("noWall", "");
    timeSelector::addOptions(true, false);

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    // Initialize name mapping table
    FieldviewNames.insert("alpha", "aalpha");
    FieldviewNames.insert("Alpha", "AAlpha");
    FieldviewNames.insert("fsmach", "ffsmach");
    FieldviewNames.insert("FSMach", "FFSMach");
    FieldviewNames.insert("re", "rre");
    FieldviewNames.insert("Re", "RRe");
    FieldviewNames.insert("time", "ttime");
    FieldviewNames.insert("Time", "TTime");
    FieldviewNames.insert("pi", "ppi");
    FieldviewNames.insert("PI", "PPI");
    FieldviewNames.insert("x", "xx");
    FieldviewNames.insert("X", "XX");
    FieldviewNames.insert("y", "yy");
    FieldviewNames.insert("Y", "YY");
    FieldviewNames.insert("z", "zz");
    FieldviewNames.insert("Z", "ZZ");
    FieldviewNames.insert("rcyl", "rrcyl");
    FieldviewNames.insert("Rcyl", "RRcyl");
    FieldviewNames.insert("theta", "ttheta");
    FieldviewNames.insert("Theta", "TTheta");
    FieldviewNames.insert("rsphere", "rrsphere");
    FieldviewNames.insert("Rsphere", "RRsphere");
    FieldviewNames.insert("k", "kk");
    FieldviewNames.insert("K", "KK");


    // Scan for all available fields, in all timesteps
    //     volScalarNames  : all scalar fields
    //     volVectorNames  : ,,  vector ,,
    //     surfScalarNames : surface fields
    //     surfVectorNames : ,,
    //     sprayScalarNames: spray fields
    //     sprayVectorNames: ,,
#   include "getFieldNames.H"

    bool hasLagrangian = false;
    if (sprayScalarNames.size() || sprayVectorNames.size())
    {
        hasLagrangian = true;
    }

    Info<< "All fields:   Foam/Fieldview" << endl;
    Info<< "    volScalar   :";
    printNames(volScalarNames, Info);
    Info<< endl;
    Info<< "    volVector   :";
    printNames(volVectorNames, Info);
    Info<< endl;
    Info<< "    surfScalar  :";
    printNames(surfScalarNames, Info);
    Info<< endl;
    Info<< "    surfVector  :";
    printNames(surfVectorNames, Info);
    Info<< endl;
    Info<< "    sprayScalar :";
    printNames(sprayScalarNames, Info);
    Info<< endl;
    Info<< "    sprayVector :";
    printNames(sprayVectorNames, Info);
    Info<< endl;


    //
    // Start writing
    //

    // make a directory called FieldView in the case
    fileName fvPath(runTime.path()/"Fieldview");

    if (isDir(fvPath))
    {
        rmDir(fvPath);
    }

    mkDir(fvPath);

    fileName fvParticleFileName(fvPath/runTime.caseName() + ".fvp");
    if (hasLagrangian)
    {
        Info<< "Opening particle file " << fvParticleFileName << endl;
    }
    std::ofstream fvParticleFile(fvParticleFileName.c_str());

    // Write spray file header
    if (hasLagrangian)
    {
#       include "writeSprayHeader.H"
    }

    // Current mesh. Start off from unloaded mesh.
    autoPtr<fieldviewTopology> topoPtr(NULL);

    label fieldViewTime = 0;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time: " << runTime.timeName() << endl;

        fvMesh::readUpdateState state = mesh.readUpdate();

        if
        (
            timeI == 0
         || state == fvMesh::TOPO_CHANGE
         || state == fvMesh::TOPO_PATCH_CHANGE
        )
        {
            // Mesh topo changed. Update Fieldview topo.

            topoPtr.reset
            (
                new fieldviewTopology
                (
                    mesh,
                    !args.optionFound("noWall")
                )
            );

            Info<< "    Mesh read:" << endl
                << "        tet   : " << topoPtr().nTet() << endl
                << "        hex   : " << topoPtr().nHex() << endl
                << "        prism : " << topoPtr().nPrism() << endl
                << "        pyr   : " << topoPtr().nPyr() << endl
                << "        poly  : " << topoPtr().nPoly() << endl
                << endl;
        }
        else if (state == fvMesh::POINTS_MOVED)
        {
            // points exists for time step, let's read them
            Info<< "    Points file detected - updating points" << endl;
        }

        const fieldviewTopology& topo = topoPtr();


        //
        // Create file and write header
        //

        fileName fvFileName
        (
            fvPath/runTime.caseName() + "_" + Foam::name(timeI) + ".uns"
        );

        Info<< "    file:" << fvFileName.c_str() << endl;


        std::ofstream fvFile(fvFileName.c_str());

        //Info<< "Writing header ..." << endl;

        // Output the magic number.
        writeInt(fvFile, FV_MAGIC);

        // Output file header and version number.
        writeStr80(fvFile, "FIELDVIEW");

        // This version of the FIELDVIEW unstructured file is "3.0".
        // This is written as two integers.
        writeInt(fvFile, 3);
        writeInt(fvFile, 0);


        // File type code. Grid and results.
        writeInt(fvFile, FV_COMBINED_FILE);

        // Reserved field, always zero
        writeInt(fvFile, 0);

        // Output constants for time, fsmach, alpha and re.
        float fBuf[4];
        fBuf[0] = runTime.value();
        fBuf[1] = 0.0;
        fBuf[2] = 0.0;
        fBuf[3] = 1.0;
        fvFile.write(reinterpret_cast<char*>(fBuf), 4*sizeof(float));


        // Output the number of grids
        writeInt(fvFile, 1);


        //
        //  Boundary table
        //
        //Info<< "Writing boundary table ..." << endl;

        // num patches
        writeInt(fvFile, mesh.boundary().size());

        forAll (mesh.boundary(), patchI)
        {
            const fvPatch& currPatch = mesh.boundary()[patchI];

            writeInt(fvFile, 1);   // data present
            writeInt(fvFile, 1);   // normals ok

            // name
            writeStr80(fvFile, currPatch.name().c_str());
        }


        //
        // Create fields:
        //     volFieldPtrs     : List of ptrs to all volScalar/Vector fields
        //                        (null if field not present at current time)
        //     volFieldNames    : FieldView compatible names of volFields
        //     surfFieldPtrs    : same for surfaceScalar/Vector
        //     surfFieldNames
#       include "createFields.H"



        //
        // Write Variables table
        //

        //Info<< "Writing variables table ..." << endl;

        writeInt(fvFile, volFieldNames.size());
        forAll(volFieldNames, fieldI)
        {
            if (volFieldPtrs[fieldI] == NULL)
            {
                Info<< "    dummy field for "
                    << volFieldNames[fieldI].c_str() << endl;
            }

            writeStr80(fvFile, volFieldNames[fieldI].c_str());
        }

        //
        // Write Boundary Variables table = vol + surface fields
        //

        //Info<< "Writing boundary variables table ..." << endl;

        writeInt
        (
            fvFile,
            volFieldNames.size() + surfFieldNames.size()
        );
        forAll(volFieldNames, fieldI)
        {
            writeStr80(fvFile, volFieldNames[fieldI].c_str());
        }
        forAll(surfFieldNames, fieldI)
        {
            if (surfFieldPtrs[fieldI] == NULL)
            {
                Info<< "    dummy surface field for "
                    << surfFieldNames[fieldI].c_str() << endl;
            }

            writeStr80(fvFile, surfFieldNames[fieldI].c_str());
        }


        // Output grid data.

        //
        // Nodes
        //

        //Info<< "Writing points ..." << endl;

        const pointField& points = mesh.points();
        label nPoints = points.size();

        writeInt(fvFile, FV_NODES);
        writeInt(fvFile, nPoints);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            floatField fField(nPoints);

            for (label pointi = 0; pointi<nPoints; pointi++)
            {
                fField[pointi] = float(points[pointi][cmpt]);
            }

            fvFile.write
            (
                reinterpret_cast<char*>(fField.begin()),
                fField.size()*sizeof(float)
            );
        }

        //
        // Boundary Faces - regular
        //

        //Info<< "Writing regular boundary faces ..." << endl;

        forAll(mesh.boundary(), patchI)
        {
            label nQuadFaces = topo.quadFaceLabels()[patchI].size()/4;

            if (nQuadFaces != 0)
            {
                writeInt(fvFile, FV_FACES);
                writeInt(fvFile, patchI + 1);  // patch number
                writeInt(fvFile, nQuadFaces);  // number of faces in patch
                fvFile.write
                (
                    reinterpret_cast<const char*>
                        (topo.quadFaceLabels()[patchI].begin()),
                    nQuadFaces*4*sizeof(int)
                );
            }
        }

        //
        // Boundary Faces - arbitrary polygon
        //

        //Info<< "Write polygonal boundary faces ..." << endl;

        forAll(mesh.boundary(), patchI)
        {
            if (topo.nPolyFaces()[patchI] > 0)
            {
                writeInt(fvFile, FV_ARB_POLY_FACES);
                writeInt(fvFile, patchI + 1);

                // number of arb faces in patch
                writeInt(fvFile, topo.nPolyFaces()[patchI]);

                const polyPatch& patchFaces = mesh.boundary()[patchI].patch();

                forAll(patchFaces, faceI)
                {
                    const face& f = patchFaces[faceI];

                    if (f.size() > 4)
                    {
                        writeInt(fvFile, f.size());

                        forAll(f, fp)
                        {
                            writeInt(fvFile, f[fp] + 1);
                        }
                    }
                }
            }
        }


        //
        // Write regular topology
        //

        //Info<< "Writing regular elements ..." << endl;

        writeInt(fvFile, FV_ELEMENTS);
        writeInt(fvFile, topo.nTet());
        writeInt(fvFile, topo.nHex());
        writeInt(fvFile, topo.nPrism());
        writeInt(fvFile, topo.nPyr());
        fvFile.write
        (
            reinterpret_cast<const char*>(topo.tetLabels().begin()),
            topo.nTet()*(1+4)*sizeof(int)
        );
        fvFile.write
        (
            reinterpret_cast<const char*>(topo.hexLabels().begin()),
            topo.nHex()*(1+8)*sizeof(int)
        );
        fvFile.write
        (
            reinterpret_cast<const char*>(topo.prismLabels().begin()),
            topo.nPrism()*(1+6)*sizeof(int)
        );
        fvFile.write
        (
            reinterpret_cast<const char*>(topo.pyrLabels().begin()),
            topo.nPyr()*(1+5)*sizeof(int)
        );


        //
        // Write arbitrary polyhedra
        //

        //Info<< "Writing polyhedral elements ..." << endl;


        const cellShapeList& cellShapes = mesh.cellShapes();
        const cellModel& unknown = *(cellModeller::lookup("unknown"));

        if (topo.nPoly() > 0)
        {
            writeInt(fvFile, FV_ARB_POLY_ELEMENTS);
            writeInt(fvFile, topo.nPoly());

            forAll(cellShapes, celli)
            {
                if (cellShapes[celli].model() == unknown)
                {
                    const cell& cll = mesh.cells()[celli];

                    // number of faces
                    writeInt(fvFile, cll.size());
                    // number of vertices used (no cell centre)
                    writeInt(fvFile, countVerts(mesh, celli));
                    // cell centre node id
                    writeInt(fvFile, -1);

                    forAll(cll, cellFacei)
                    {
                        label faceI = cll[cellFacei];

                        const face& f = mesh.faces()[faceI];

                        // Not a wall for now
                        writeInt(fvFile, NOT_A_WALL);

                        writeInt(fvFile, f.size());

                        if (mesh.faceOwner()[faceI] == celli)
                        {
                            forAll(f, fp)
                            {
                                writeInt(fvFile, f[fp]+1);
                            }
                        }
                        else
                        {
                            for (label fp = f.size()-1; fp >= 0; fp--)
                            {
                                writeInt(fvFile, f[fp]+1);
                            }
                        }

                        // Number of hanging nodes
                        writeInt(fvFile, 0);
                    }
                }
            }
        }


        //
        // Variables data
        //

        //Info<< "Writing variables data ..." << endl;

        volPointInterpolation pInterp(mesh);

        writeInt(fvFile, FV_VARIABLES);


        forAll(volFieldPtrs, fieldI)
        {
            if (volFieldPtrs[fieldI] != NULL)
            {
                const volScalarField& vField = *volFieldPtrs[fieldI];

                // Interpolate to points
                pointScalarField psf(pInterp.interpolate(vField));

                floatField fField(nPoints);

                for (label pointi = 0; pointi<nPoints; pointi++)
                {
                    fField[pointi] = float(psf[pointi]);
                }

                fvFile.write
                (
                    reinterpret_cast<char*>(fField.begin()),
                    fField.size()*sizeof(float)
                );
            }
            else
            {
                // Create dummy field
                floatField dummyField(nPoints, 0.0);

                fvFile.write
                (
                    reinterpret_cast<char*>(dummyField.begin()),
                    dummyField.size()*sizeof(float)
                );
            }
        }


        //
        // Boundary variables data
        //     1. volFields
        //     2. surfFields

        //Info<< "Writing regular boundary elements data ..." << endl;

        writeInt(fvFile, FV_BNDRY_VARS);

        forAll(volFieldPtrs, fieldI)
        {
            if (volFieldPtrs[fieldI] != NULL)
            {
                const volScalarField& vsf = *volFieldPtrs[fieldI];

                forAll (mesh.boundary(), patchI)
                {
                    writeFaceData
                    (
                        mesh,
                        topo,
                        patchI,
                        vsf.boundaryField()[patchI],
                        false,
                        fvFile
                    );
                }
            }
            else
            {
                forAll (mesh.boundaryMesh(), patchI)
                {
                    // Dummy value.
                    floatField fField
                    (
                        mesh.boundaryMesh()[patchI].size()
                      - topo.nPolyFaces()[patchI],
                        0.0
                    );

                    fvFile.write
                    (
                        reinterpret_cast<char*>(fField.begin()),
                        fField.size()*sizeof(float)
                    );
                }
            }
        }

        // surfFields
        forAll(surfFieldPtrs, fieldI)
        {
            if (surfFieldPtrs[fieldI] != NULL)
            {
                const surfaceScalarField& ssf = *surfFieldPtrs[fieldI];

                forAll (mesh.boundary(), patchI)
                {
                    writeFaceData
                    (
                        mesh,
                        topo,
                        patchI,
                        ssf.boundaryField()[patchI],
                        false,
                        fvFile
                    );
                }
            }
            else
            {
                forAll (mesh.boundaryMesh(), patchI)
                {
                    // Dummy value.
                    floatField fField
                    (
                        mesh.boundaryMesh()[patchI].size()
                      - topo.nPolyFaces()[patchI],
                        0.0
                    );

                    fvFile.write
                    (
                        reinterpret_cast<char*>(fField.begin()),
                        fField.size()*sizeof(float)
                    );
                }
            }
        }

        //
        // Polygonal faces boundary data
        //     1. volFields
        //     2. surfFields

        //Info<< "Writing polygonal boundary elements data ..." << endl;

        writeInt(fvFile, FV_ARB_POLY_BNDRY_VARS);
        forAll(volFieldPtrs, fieldI)
        {
            if (volFieldPtrs[fieldI] != NULL)
            {
                const volScalarField& vsf = *volFieldPtrs[fieldI];

                // All non-empty patches
                forAll (mesh.boundary(), patchI)
                {
                    writeFaceData
                    (
                        mesh,
                        topo,
                        patchI,
                        vsf.boundaryField()[patchI],
                        true,
                        fvFile
                    );
                }
            }
            else
            {
                forAll (mesh.boundary(), patchI)
                {
                    // Dummy value.
                    floatField fField(topo.nPolyFaces()[patchI], 0.0);

                    fvFile.write
                    (
                        reinterpret_cast<char*>(fField.begin()),
                        fField.size()*sizeof(float)
                    );
                }
            }
        }

        // surfFields
        forAll(surfFieldPtrs, fieldI)
        {
            if (surfFieldPtrs[fieldI] != NULL)
            {
                const surfaceScalarField& ssf = *surfFieldPtrs[fieldI];

                // All non-empty patches
                forAll(mesh.boundary(), patchI)
                {
                    writeFaceData
                    (
                        mesh,
                        topo,
                        patchI,
                        ssf.boundaryField()[patchI],
                        true,
                        fvFile
                    );
                }
            }
            else
            {
                forAll (mesh.boundaryMesh(), patchI)
                {
                    // Dummy value.
                    floatField fField
                    (
                        mesh.boundaryMesh()[patchI].size()
                      - topo.nPolyFaces()[patchI],
                        0.0
                    );

                    fvFile.write
                    (
                        reinterpret_cast<char*>(fField.begin()),
                        fField.size()*sizeof(float)
                    );
                }
            }
        }


        //
        // Cleanup volume and surface fields
        //
        forAll(volFieldPtrs, fieldI)
        {
            delete volFieldPtrs[fieldI];
        }
        forAll(surfFieldPtrs, fieldI)
        {
            delete surfFieldPtrs[fieldI];
        }




        //
        // Spray
        //
        if (hasLagrangian)
        {
            // Read/create fields:
            //     sprayScalarFieldPtrs: List of ptrs to lagrangian scalfields
            //     sprayVectorFieldPtrs:               ,,           vec  ,,
#           include "createSprayFields.H"


            // Write time header

            // Time index (FieldView: has to start from 1)
            writeInt(fvParticleFile, fieldViewTime + 1);

            // Time value
            writeFloat(fvParticleFile, runTime.value());

            // Read particles
            Cloud<passiveParticle> parcels(mesh);

            // Num particles
            writeInt(fvParticleFile, parcels.size());

            Info<< "    Writing " << parcels.size() << " particles." << endl;


            //
            // Output data parcelwise
            //

            label parcelNo = 0;


            for
            (
                Cloud<passiveParticle>::iterator elmnt = parcels.begin();
                elmnt != parcels.end();
                ++elmnt, parcelNo++
            )
            {
                writeInt(fvParticleFile, parcelNo+1);

                writeFloat(fvParticleFile, elmnt().position().x());
                writeFloat(fvParticleFile, elmnt().position().y());
                writeFloat(fvParticleFile, elmnt().position().z());

                forAll(sprayScalarFieldPtrs, fieldI)
                {
                    if (sprayScalarFieldPtrs[fieldI] != NULL)
                    {
                        const IOField<scalar>& sprayField =
                            *sprayScalarFieldPtrs[fieldI];
                        writeFloat
                        (
                            fvParticleFile,
                            sprayField[parcelNo]
                        );
                    }
                    else
                    {
                        writeFloat(fvParticleFile, 0.0);
                    }
                }
                forAll(sprayVectorFieldPtrs, fieldI)
                {
                    if (sprayVectorFieldPtrs[fieldI] != NULL)
                    {
                        const IOField<vector>& sprayVectorField =
                            *sprayVectorFieldPtrs[fieldI];
                        const vector& val =
                            sprayVectorField[parcelNo];

                        writeFloat(fvParticleFile, val.x());
                        writeFloat(fvParticleFile, val.y());
                        writeFloat(fvParticleFile, val.z());
                    }
                    else
                    {
                        writeFloat(fvParticleFile, 0.0);
                        writeFloat(fvParticleFile, 0.0);
                        writeFloat(fvParticleFile, 0.0);
                    }
                }
            }

            // increment fieldView particle time
            fieldViewTime++;


            //
            // Cleanup spray fields
            //
            forAll(sprayScalarFieldPtrs, fieldI)
            {
                delete sprayScalarFieldPtrs[fieldI];
            }
            forAll(sprayVectorFieldPtrs, fieldI)
            {
                delete sprayVectorFieldPtrs[fieldI];
            }

        } // end of hasLagrangian
    }

    if (!hasLagrangian)
    {
        rm(fvParticleFileName);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
