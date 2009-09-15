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
    Reader module for Fieldview9 to read OpenFOAM mesh and data.

    Creates new 'fvbin' type executable which needs to be installed in place
    of bin/fvbin.

    Implements a reader for combined mesh&results on an unstructured mesh.

    See Fieldview Release 9 Reference Manual and coding in user/ directory
    of the Fieldview release.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOobjectList.H"
#include "GeometricField.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "readerDatabase.H"
#include "wallPolyPatch.H"
#include "ListOps.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern "C"
{
// Define various Fieldview constants and prototypes

#include "fv_reader_tags.h"

static const int FVHEX     = 2;
static const int FVPRISM   = 3;
static const int FVPYRAMID = 4;
static const int FVTET     = 1;
static const int ITYPE     = 1;

unsigned int fv_encode_elem_header(int elem_type, int wall_info[]);
void time_step_get_value(float*, int, int*, float*, int*);
void fv_error_msg(const char*, const char*);

void reg_single_unstruct_reader
(
    char *,
    void
    (
        char*, int*, int*, int*, int*, int*, int*,
        int[], int*, char[][80], int[], int*, char[][80], int*
    ),
    void
    (
        int*, int*, int*, float[], int*, float[], int*
    )
);

int create_tet(const int, const int[8], const int[]);
int create_pyramid(const int, const int[8], const int[]);
int create_prism(const int, const int[8], const int[]);
int create_hex(const int, const int[8], const int[]);

typedef unsigned char uChar;
extern uChar create_bndry_face_buffered
(
    int bndry_type,
    int num_verts,
    int verts[],
    int *normals_flags,
    int num_grid_nodes
);

/*
 * just define empty readers here for ease of linking.
 * Comment out if you have doubly defined linking error on this symbol
 */
void ftn_register_data_readers()
{}

/*
 * just define empty readers here for ease of linking.
 * Comment out if you have doubly defined linking error on this symbol
 */
void ftn_register_functions()
{}

/*
 * just define empty readers here for ease of linking.
 * Comment out if you have doubly defined linking error on this symbol
 */
void register_functions()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//
// Storage for all Foam state (mainly database & mesh)
//
static readerDatabase db_;


// Write fv error message.
static void errorMsg(const string& msg)
{
    fv_error_msg("Foam Reader", msg.c_str());
}


// Simple check if directory is valid case directory.
static bool validCase(const fileName& rootAndCase)
{
    //if (isDir(rootAndCase/"system") && isDir(rootAndCase/"constant"))
    if (isDir(rootAndCase/"constant"))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Check whether case has topo changes by looking back from last time
// to first directory with polyMesh/cells.
static bool hasTopoChange(const instantList& times)
{
    label lastIndex = times.size()-1;

    const Time& runTime = db_.runTime();

    // Only set time; do not update mesh.
    runTime.setTime(times[lastIndex], lastIndex);

    fileName facesInst(runTime.findInstance(polyMesh::meshSubDir, "faces"));

    // See if cellInst is constant directory. Note extra .name() is for
    // parallel cases when runTime.constant is '../constant'
    if (facesInst != runTime.constant().name())
    {
        Info<< "Found cells file in " << facesInst << " so case has "
            << "topological changes" << endl;

        return true;
    }
    else
    {
        Info<< "Found cells file in " << facesInst << " so case has "
            << "no topological changes" << endl;

        return false;
    }
}


static bool selectTime(const instantList& times, int* iret)
{
    List<float> fvTimes(2*times.size());

    forAll(times, timeI)
    {
        fvTimes[2*timeI]   = float(timeI);
        fvTimes[2*timeI+1] = float(times[timeI].value());
    }

    int istep;

    float time;

    *iret=0;

    time_step_get_value(fvTimes.begin(), times.size(), &istep, &time, iret);

    if (*iret == -5)
    {
        errorMsg("Out of memory.");

        return false;
    }
    if (*iret == -15)
    {
        // Cancel action.
        return false;
    }
    if (*iret != 0)
    {
        errorMsg("Unspecified error.");

        return false;
    }
    Info<< "Selected timeStep:" << istep << "  time:" << scalar(time) << endl;

    // Set time and load mesh.
    db_.setTime(times[istep], istep);

    return true;
}


// Gets (names of) all fields in all timesteps.
static void createFieldNames
(
    const Time& runTime,
    const instantList& Times,
    const word& setName
)
{
    // From foamToFieldView9/getFieldNames.H:

    HashSet<word> volScalarHash;
    HashSet<word> volVectorHash;
    HashSet<word> surfScalarHash;
    HashSet<word> surfVectorHash;

    if (setName.empty())
    {
        forAll(Times, timeI)
        {
            const word& timeName = Times[timeI].name();

            // Add all fields to hashtable
            IOobjectList objects(runTime, timeName);

            wordList vsNames(objects.names(volScalarField::typeName));

            forAll(vsNames, fieldI)
            {
                volScalarHash.insert(vsNames[fieldI]);
            }

            wordList vvNames(objects.names(volVectorField::typeName));

            forAll(vvNames, fieldI)
            {
                volVectorHash.insert(vvNames[fieldI]);
            }
        }
    }
    db_.setFieldNames(volScalarHash.toc(), volVectorHash.toc());
}


// Appends interpolated values of fieldName to vars array.
static void storeScalarField
(
    const volPointInterpolation& pInterp,
    const word& fieldName,
    float vars[],
    label& pointI
)
{
    label nPoints = db_.mesh().nPoints();
    label nTotPoints = nPoints + db_.polys().size();

    // Check if present
    IOobject ioHeader
    (
        fieldName,
        db_.runTime().timeName(),
        db_.runTime(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (ioHeader.headerOk())
    {
        Info<< "Storing " << nTotPoints << " of interpolated " << fieldName
            << endl;

        volScalarField field(ioHeader, db_.mesh());

        pointScalarField psf(pInterp.interpolate(field));

        forAll(psf, i)
        {
            vars[pointI++] = float(psf[i]);
        }

        const labelList& polys = db_.polys();

        forAll(polys, i)
        {
            label cellI = polys[i];

            vars[pointI++] = float(field[cellI]);
        }
    }
    else
    {
        Info<< "Storing " << nTotPoints << " of dummy values of " << fieldName
            << endl;

        for(label i = 0; i < nPoints; i++)
        {
            vars[pointI++] = 0.0;
        }

        const labelList& polys = db_.polys();

        forAll(polys, i)
        {
            vars[pointI++] = 0.0;
        }
    }
}


// Appends interpolated values of fieldName to vars array.
static void storeVectorField
(
    const volPointInterpolation& pInterp,
    const word& fieldName,
    float vars[],
    label& pointI
)
{
    label nPoints = db_.mesh().nPoints();

    label nTotPoints = nPoints + db_.polys().size();

    // Check if present
    IOobject ioHeader
    (
        fieldName,
        db_.runTime().timeName(),
        db_.runTime(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (ioHeader.headerOk())
    {
        Info<< "Storing " << nTotPoints << " of interpolated " << fieldName
            << endl;

        volVectorField field(ioHeader, db_.mesh());

        for (direction d = 0; d < vector::nComponents; d++)
        {
            tmp<volScalarField> tcomp = field.component(d);
            const volScalarField& comp = tcomp();

            pointScalarField psf(pInterp.interpolate(comp));

            forAll(psf, i)
            {
                vars[pointI++] = float(psf[i]);
            }

            const labelList& polys = db_.polys();

            forAll(polys, i)
            {
                label cellI = polys[i];

                vars[pointI++] = float(comp[cellI]);
            }
        }
    }
    else
    {
        Info<< "Storing " << nTotPoints << " of dummy values of " << fieldName
            << endl;

        for (direction d = 0; d < vector::nComponents; d++)
        {
            for(label i = 0; i < nPoints; i++)
            {
                vars[pointI++] = 0.0;
            }

            const labelList& polys = db_.polys();

            forAll(polys, i)
            {
                vars[pointI++] = 0.0;
            }
        }
    }
}


// Returns Fieldview face_type of mesh face faceI.
static label getFvType(const polyMesh& mesh, const label faceI)
{
    return mesh.boundaryMesh().whichPatch(faceI) + 1;
}


// Returns Fieldview face_type of face f.
static label getFaceType
(
    const polyMesh& mesh,
    const labelList& faceLabels,
    const face& f
)
{
    // Search in subset faceLabels of faces for index of face f.
    const faceList& faces = mesh.faces();

    forAll(faceLabels, i)
    {
        label faceI = faceLabels[i];

        if (f == faces[faceI])
        {
            // Convert patch to Fieldview face_type.
            return getFvType(mesh, faceI);
        }
    }

    FatalErrorIn("getFaceType")
        << "Cannot find face " << f << " in mesh face subset " << faceLabels
        << abort(FatalError);

    return -1;
}


// Returns Fieldview face_types for set of faces
static labelList getFaceTypes
(
    const polyMesh& mesh,
    const labelList& cellFaces,
    const faceList& cellShapeFaces
)
{
    labelList faceLabels(cellShapeFaces.size());

    forAll(cellShapeFaces, i)
    {
        faceLabels[i] = getFaceType(mesh, cellFaces, cellShapeFaces[i]);
    }
    return faceLabels;
}


/*
 * Callback for querying file contents. Taken from user/user_unstruct_combined.f
 */
void user_query_file_function
(
    /* input */
    char* fname,        /* filename */
    int* lenf,          /* length of fName */
    int* iunit,         /* fortran unit to use */
    int* max_grids,     /* maximum number of grids allowed */
    int* max_face_types,/* maximum number of face types allowed */
    int* max_vars,      /* maximum number of result variables allowed per */
                        /* grid point*/

    /* output */
    int* num_grids,     /* number of grids that will be read from data file */
    int  num_nodes[],   /* (array of node counts for all grids) */
    int* num_face_types,        /* number of special face types */
    char face_type_names[][80], /* array of face-type names */
    int  wall_flags[],          /* array of flags for the "wall" behavior */
    int* num_vars,              /* number of result variables per grid point */
    char var_names[][80],       /* array of variable names */
    int* iret                   /* return value */
)
{
    fprintf(stderr, "\n** user_query_file_function\n");

    string rootAndCaseString(fname, *lenf);

    fileName rootAndCase(rootAndCaseString);

    word setName("");

    if (!validCase(rootAndCase))
    {
        setName = rootAndCase.name();

        rootAndCase = rootAndCase.path();

        word setDir = rootAndCase.name();

        rootAndCase = rootAndCase.path();

        word meshDir = rootAndCase.name();

        rootAndCase = rootAndCase.path();
        rootAndCase = rootAndCase.path();

        if
        (
            setDir == "sets"
         && meshDir == polyMesh::typeName
         && validCase(rootAndCase)
        )
        {
            // Valid set (hopefully - cannot check contents of setName yet).
        }
        else
        {
            errorMsg
            (
                "Could not find system/ and constant/ directory in\n"
              + rootAndCase
              + "\nPlease select a Foam case directory."
            );

            *iret = 1;

            return;
        }

    }

    fileName rootDir(rootAndCase.path());

    fileName caseName(rootAndCase.name());

    // handle trailing '/'
    if (caseName.empty())
    {
        caseName = rootDir.name();
        rootDir  = rootDir.path();
    }

    Info<< "rootDir  : " << rootDir << endl
        << "caseName : " << caseName << endl
        << "setName  : " << setName << endl;

    //
    // Get/reuse database and mesh
    //

    bool caseChanged = db_.setRunTime(rootDir, caseName, setName);


    //
    // Select time
    //

    instantList Times = db_.runTime().times();

    // If topo case set database time and update mesh.
    if (hasTopoChange(Times))
    {
        if (!selectTime(Times, iret))
        {
            return;
        }
    }
    else if (caseChanged)
    {
        // Load mesh (if case changed) to make sure we have nPoints etc.
        db_.loadMesh();
    }


    //
    // Set output variables
    //

    *num_grids = 1;

    const fvMesh& mesh = db_.mesh();

    label nTotPoints = mesh.nPoints() + db_.polys().size();

    num_nodes[0] = nTotPoints;

    Info<< "setting num_nodes:" << num_nodes[0] << endl;

    Info<< "setting num_face_types:" << mesh.boundary().size() << endl;

    *num_face_types = mesh.boundary().size();

    if (*num_face_types > *max_face_types)
    {
        errorMsg("Too many patches. FieldView limit:" + name(*max_face_types));

        *iret = 1;

        return;
    }


    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        strcpy(face_type_names[patchI], patch.name().c_str());

        if (isA<wallPolyPatch>(patch))
        {
            wall_flags[patchI] = 1;
        }
        else
        {
            wall_flags[patchI] = 0;
        }
        Info<< "Patch " << patch.name() << " is wall:"
            <<  wall_flags[patchI] << endl;
    }

    //- Find all volFields and add them to database
    createFieldNames(db_.runTime(), Times, setName);

    *num_vars = db_.volScalarNames().size() + 3*db_.volVectorNames().size();

    if (*num_vars > *max_vars)
    {
        errorMsg("Too many variables. FieldView limit:" + name(*max_vars));

        *iret = 1;

        return;
    }


    int nameI = 0;

    forAll(db_.volScalarNames(), i)
    {
        const word& fieldName = db_.volScalarNames()[i];

        const word& fvName = db_.getFvName(fieldName);

        strcpy(var_names[nameI++], fvName.c_str());
    }

    forAll(db_.volVectorNames(), i)
    {
        const word& fieldName = db_.volVectorNames()[i];

        const word& fvName = db_.getFvName(fieldName);

        strcpy(var_names[nameI++], (fvName + "x;" + fvName).c_str());
        strcpy(var_names[nameI++], (fvName + "y").c_str());
        strcpy(var_names[nameI++], (fvName + "z").c_str());
    }

    *iret = 0;
}


/*
 * Callback for reading timestep. Taken from user/user_unstruct_combined.f
 */
void user_read_one_grid_function
(
    int* iunit,     /* in: fortran unit number */
    int* igrid,     /* in: grid number to read */
    int* nodecnt,   /* in: number of nodes to read */
    float xyz[],    /* out: coordinates of nodes: x1..xN y1..yN z1..zN */
    int* num_vars,  /* in: number of results per node */
    float vars[],   /* out: values per node */
    int* iret       /* out: return value */
)
{
    fprintf(stderr, "\n** user_read_one_grid_function\n");

    if (*igrid != 1)
    {
        errorMsg("Illegal grid number " + Foam::name(*igrid));

        *iret = 1;

        return;
    }

    // Get current time
    instantList Times = db_.runTime().times();

    // Set database time and update mesh.
    // Note: this should not be nessecary here. We already have the correct
    // time set and mesh loaded. This is only nessecary because Fieldview
    // otherwise thinks the case is non-transient.
    if (!selectTime(Times, iret))
    {
        return;
    }


    const fvMesh& mesh = db_.mesh();

    // With mesh now loaded check for change in number of points.
    label nTotPoints = mesh.nPoints() + db_.polys().size();

    if (*nodecnt != nTotPoints)
    {
        errorMsg
        (
            "nodecnt differs from number of points in mesh.\nnodecnt:"
          + Foam::name(*nodecnt)
          + "  mesh:"
          + Foam::name(nTotPoints)
        );

        *iret = 1;

        return;
    }


    if
    (
        *num_vars
     != (db_.volScalarNames().size() + 3*db_.volVectorNames().size())
    )
    {
        errorMsg("Illegal number of variables " + name(*num_vars));

        *iret = 1;

        return;
    }

    //
    // Set coordinates
    //

    const pointField& points = mesh.points();

    int xIndex = 0;
    int yIndex = xIndex + nTotPoints;
    int zIndex = yIndex + nTotPoints;

    // Add mesh points first.
    forAll(points, pointI)
    {
        xyz[xIndex++] = points[pointI].x();
        xyz[yIndex++] = points[pointI].y();
        xyz[zIndex++] = points[pointI].z();
    }

    // Add cell centres of polys
    const pointField& ctrs = mesh.cellCentres();

    const labelList& polys = db_.polys();

    forAll(polys, i)
    {
        label cellI = polys[i];

        xyz[xIndex++] = ctrs[cellI].x();
        xyz[yIndex++] = ctrs[cellI].y();
        xyz[zIndex++] = ctrs[cellI].z();
    }


    //
    // Define elements by calling fv routines
    //

    static const cellModel& tet = *(cellModeller::lookup("tet"));
    static const cellModel& tetWedge = *(cellModeller::lookup("tetWedge"));
    static const cellModel& pyr = *(cellModeller::lookup("pyr"));
    static const cellModel& prism = *(cellModeller::lookup("prism"));
    static const cellModel& wedge = *(cellModeller::lookup("wedge"));
    static const cellModel& hex = *(cellModeller::lookup("hex"));
    //static const cellModel& splitHex = *(cellModeller::lookup("splitHex"));

    int tetVerts[4];
    int pyrVerts[5];
    int prismVerts[6];
    int hexVerts[8];

    int tetFaces[4];
    int pyrFaces[5];
    int prismFaces[5];
    int hexFaces[6];

    const cellShapeList& cellShapes = mesh.cellShapes();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelList& owner = mesh.faceOwner();
    label nPoints = mesh.nPoints();

    // Fieldview face_types array with all faces marked as internal.
    labelList internalFaces(6, 0);

    // Mark all cells next to boundary so we don't have to calculate
    // wall_types for internal cells and can just pass internalFaces.
    boolList wallCell(mesh.nCells(), false);

    label nWallCells = 0;

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        label cellI = owner[faceI];

        if (!wallCell[cellI])
        {
            wallCell[cellI] = true;

            nWallCells++;
        }
    }

    label nPolys = 0;

    forAll(cellShapes, cellI)
    {
        const cellShape& cellShape = cellShapes[cellI];
        const cellModel& cellModel = cellShape.model();
        const cell& cellFaces = cells[cellI];

        int istat = 0;

        if (cellModel == tet)
        {
            tetVerts[0] = cellShape[3] + 1;
            tetVerts[1] = cellShape[0] + 1;
            tetVerts[2] = cellShape[1] + 1;
            tetVerts[3] = cellShape[2] + 1;

            if (wallCell[cellI])
            {
                labelList faceTypes =
                    getFaceTypes(mesh, cellFaces, cellShape.faces());

                tetFaces[0] = faceTypes[2];
                tetFaces[1] = faceTypes[3];
                tetFaces[2] = faceTypes[0];
                tetFaces[3] = faceTypes[1];

                istat = create_tet(ITYPE, tetVerts, tetFaces);
            }
            else
            {
                // All faces internal so use precalculated zero.
                istat = create_tet(ITYPE, tetVerts, internalFaces.begin());
            }
        }
        else if (cellModel == tetWedge)
        {
            prismVerts[0] = cellShape[0] + 1;
            prismVerts[1] = cellShape[3] + 1;
            prismVerts[2] = cellShape[4] + 1;
            prismVerts[3] = cellShape[1] + 1;
            prismVerts[4] = cellShape[4] + 1;
            prismVerts[5] = cellShape[2] + 1;

            if (wallCell[cellI])
            {
                labelList faceTypes =
                    getFaceTypes(mesh, cellFaces, cellShape.faces());

                prismFaces[0] = faceTypes[1];
                prismFaces[1] = faceTypes[2];
                prismFaces[2] = faceTypes[3];
                prismFaces[3] = faceTypes[0];
                prismFaces[4] = faceTypes[3];

                istat = create_prism(ITYPE, prismVerts, prismFaces);
            }
            else
            {
                istat = create_prism(ITYPE, prismVerts, internalFaces.begin());
            }
        }
        else if (cellModel == pyr)
        {
            pyrVerts[0] = cellShape[0] + 1;
            pyrVerts[1] = cellShape[1] + 1;
            pyrVerts[2] = cellShape[2] + 1;
            pyrVerts[3] = cellShape[3] + 1;
            pyrVerts[4] = cellShape[4] + 1;

            if (wallCell[cellI])
            {
                labelList faceTypes =
                    getFaceTypes(mesh, cellFaces, cellShape.faces());

                pyrFaces[0] = faceTypes[0];
                pyrFaces[1] = faceTypes[3];
                pyrFaces[2] = faceTypes[2];
                pyrFaces[3] = faceTypes[1];
                pyrFaces[4] = faceTypes[4];

                istat = create_pyramid(ITYPE, pyrVerts, pyrFaces);
            }
            else
            {
                istat = create_pyramid(ITYPE, pyrVerts, internalFaces.begin());
            }
        }
        else if (cellModel == prism)
        {
            prismVerts[0] = cellShape[0] + 1;
            prismVerts[1] = cellShape[3] + 1;
            prismVerts[2] = cellShape[4] + 1;
            prismVerts[3] = cellShape[1] + 1;
            prismVerts[4] = cellShape[5] + 1;
            prismVerts[5] = cellShape[2] + 1;

            if (wallCell[cellI])
            {
                labelList faceTypes =
                    getFaceTypes(mesh, cellFaces, cellShape.faces());

                prismFaces[0] = faceTypes[4];
                prismFaces[1] = faceTypes[2];
                prismFaces[2] = faceTypes[3];
                prismFaces[3] = faceTypes[0];
                prismFaces[4] = faceTypes[1];

                istat = create_prism(ITYPE, prismVerts, prismFaces);
            }
            else
            {
                istat = create_prism(ITYPE, prismVerts, internalFaces.begin());
            }
        }
        else if (cellModel == wedge)
        {
            hexVerts[0] = cellShape[0] + 1;
            hexVerts[1] = cellShape[1] + 1;
            hexVerts[2] = cellShape[0] + 1;
            hexVerts[3] = cellShape[2] + 1;
            hexVerts[4] = cellShape[3] + 1;
            hexVerts[5] = cellShape[4] + 1;
            hexVerts[6] = cellShape[6] + 1;
            hexVerts[7] = cellShape[5] + 1;

            if (wallCell[cellI])
            {
                labelList faceTypes =
                    getFaceTypes(mesh, cellFaces, cellShape.faces());

                hexFaces[0] = faceTypes[2];
                hexFaces[1] = faceTypes[3];
                hexFaces[2] = faceTypes[0];
                hexFaces[3] = faceTypes[1];
                hexFaces[4] = faceTypes[4];
                hexFaces[5] = faceTypes[5];

                istat = create_hex(ITYPE, hexVerts, hexFaces);
            }
            else
            {
                istat = create_hex(ITYPE, hexVerts, internalFaces.begin());
            }
        }
        else if (cellModel == hex)
        {
            hexVerts[0] = cellShape[0] + 1;
            hexVerts[1] = cellShape[1] + 1;
            hexVerts[2] = cellShape[3] + 1;
            hexVerts[3] = cellShape[2] + 1;
            hexVerts[4] = cellShape[4] + 1;
            hexVerts[5] = cellShape[5] + 1;
            hexVerts[6] = cellShape[7] + 1;
            hexVerts[7] = cellShape[6] + 1;

            if (wallCell[cellI])
            {
                labelList faceTypes =
                    getFaceTypes(mesh, cellFaces, cellShape.faces());

                hexFaces[0] = faceTypes[0];
                hexFaces[1] = faceTypes[1];
                hexFaces[2] = faceTypes[4];
                hexFaces[3] = faceTypes[5];
                hexFaces[4] = faceTypes[2];
                hexFaces[5] = faceTypes[3];

                istat = create_hex(ITYPE, hexVerts, hexFaces);
            }
            else
            {
                istat = create_hex(ITYPE, hexVerts, internalFaces.begin());
            }
        }
        else
        {
            forAll(cellFaces, cFaceI)
            {
                label faceI = cellFaces[cFaceI];

                // Get Fieldview facetype (internal/on patch)
                label fvFaceType = getFvType(mesh, faceI);

                const face& f = faces[faceI];

                // Indices into storage
                label nQuads = 0;
                label nTris = 0;

                // Storage for triangles and quads created by face
                // decomposition (sized for worst case)
                faceList quadFaces((f.size() - 2)/2);
                faceList triFaces(f.size() - 2);

                f.trianglesQuads
                (
                    points,
                    nTris,
                    nQuads,
                    triFaces,
                    quadFaces
                );

                // Label of cell centre in fv point list.
                label polyCentrePoint = nPoints + nPolys;

                for (label i=0; i<nTris; i++)
                {
                    if (cellI == owner[faceI])
                    {
                        tetVerts[0] = triFaces[i][0] + 1;
                        tetVerts[1] = triFaces[i][1] + 1;
                        tetVerts[2] = triFaces[i][2] + 1;
                        tetVerts[3] = polyCentrePoint + 1;
                    }
                    else
                    {
                        tetVerts[0] = triFaces[i][2] + 1;
                        tetVerts[1] = triFaces[i][1] + 1;
                        tetVerts[2] = triFaces[i][0] + 1;
                        tetVerts[3] = polyCentrePoint + 1;
                    }

                    if (wallCell[cellI])
                    {
                        // Outside face is one without polyCentrePoint
                        tetFaces[0] = fvFaceType;
                        tetFaces[1] = 0;
                        tetFaces[2] = 0;
                        tetFaces[3] = 0;

                        istat = create_tet(ITYPE, tetVerts, tetFaces);
                    }
                    else
                    {
                        istat =
                            create_tet
                            (
                                ITYPE,
                                tetVerts,
                                internalFaces.begin()
                            );
                    }
                }

                for (label i=0; i<nQuads; i++)
                {
                    if (cellI == owner[faceI])
                    {
                        pyrVerts[0] = quadFaces[i][3] + 1;
                        pyrVerts[1] = quadFaces[i][2] + 1;
                        pyrVerts[2] = quadFaces[i][1] + 1;
                        pyrVerts[3] = quadFaces[i][0] + 1;
                        pyrVerts[4] = polyCentrePoint + 1;

                    }
                    else
                    {
                        pyrVerts[0] = quadFaces[i][0] + 1;
                        pyrVerts[1] = quadFaces[i][1] + 1;
                        pyrVerts[2] = quadFaces[i][2] + 1;
                        pyrVerts[3] = quadFaces[i][3] + 1;
                        pyrVerts[4] = polyCentrePoint + 1;
                    }

                    if (wallCell[cellI])
                    {
                        // Outside face is one without polyCentrePoint
                        pyrFaces[0] = fvFaceType;
                        pyrFaces[1] = 0;
                        pyrFaces[2] = 0;
                        pyrFaces[3] = 0;
                        pyrFaces[4] = 0;

                        istat = create_pyramid(ITYPE, pyrVerts, pyrFaces);
                    }
                    else
                    {
                        istat =
                            create_pyramid
                            (
                                ITYPE,
                                pyrVerts,
                                internalFaces.begin()
                            );
                    }
                }

                if (istat != 0)
                {
                    errorMsg("Error during adding cell " + name(cellI));

                    *iret = 1;

                    return;
                }
            }

            nPolys++;
        }

        if (istat != 0)
        {
            errorMsg("Error during adding cell " + name(cellI));

            *iret = 1;

            return;
        }
    }



    //
    // Set fieldvalues
    //

    pointMesh pMesh(mesh);

    volPointInterpolation pInterp(mesh, pMesh);

    int pointI = 0;

    forAll(db_.volScalarNames(), i)
    {
        const word& fieldName = db_.volScalarNames()[i];

        storeScalarField(pInterp, fieldName, vars, pointI);
    }


    forAll(db_.volVectorNames(), i)
    {
        const word& fieldName = db_.volVectorNames()[i];

        storeVectorField(pInterp, fieldName, vars, pointI);
    }

    // Return without failure.
    *iret = 0;
}


void register_data_readers()
{
    /*
    **
    ** You should edit this file to "register" a user-defined data
    ** reader with FIELDVIEW, if the user functions being registered
    ** here are written in C.
    ** You should edit "ftn_register_data_readers.f" if the user functions
    ** being registered are written in Fortran.
    ** In either case, the user functions being registered may call other
    ** functions written in either language (C or Fortran); only the
    ** language of the "query" and "read" functions referenced here matters
    ** to FIELDVIEW.
    **
    ** The following shows a sample user-defined data reader being
    ** "registered" with FIELDVIEW.
    **
    ** The "extern void" declarations should match the names of the
    ** query and grid-reading functions you are providing. It is
    ** strongly suggested that all such names begin with "user" so
    ** as not to conflict with global names in FIELDVIEW.
    **
    ** You may call any combination of the data reader registration
    ** functions shown below ("register_two_file_reader" and/or
    ** "register_single_file_reader" and/or "register_single_unstruct_reader"
    ** and/or "register_double_unstruct_reader") as many times as you like,
    ** in order to create several different data readers. Each data reader
    ** should of course have different "query" and "read" functions, all of
    ** which should also appear in "extern void" declarations.
    **
    */

    /*
    ** like this for combined unstructured grids & results in a single file
    */
    reg_single_unstruct_reader (
	"Foam Reader",		       /* title you want for data reader */
	user_query_file_function,      /* whatever you called this */
	user_read_one_grid_function    /* whatever you called this */
	);
}




}


// ************************************************************************* //
