/*--------------------------------*- C++ -*----------------------------------*\
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

Application
    ansysToFoam

Description
    Converts an ANSYS input mesh file, exported from I-DEAS, to FOAM format.

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

 /* ------------------------------------------------------------------------- *\
   ------ local definitions
 \* ------------------------------------------------------------------------- */

#include <sstream>

#include "scalar.H"
#include "IStringStream.H"

using namespace Foam;

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "preservePatchTypes.H"
#include "cellShape.H"
#include "cellModeller.H"
#include "SLList.H"
#include "SLPtrList.H"

SLList<point> slPoints;
SLList<label> slPointMap;
label maxNodei = 0;

SLPtrList<labelList> slCellLabels;
SLList<label> slCellMap;
label maxCelli = 0;

PtrList<SLList<label> > slPatchCells;
PtrList<SLList<label> > slPatchCellFaces;


// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
#if YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+
cspace                ","{space}

alpha                 [_A-Za-z]
digit                 [0-9]
dec_digit             [0-9]
octal_digit           [0-7]
hex_digit             [0-9a-fA-F]

identifier            {alpha}({alpha}|{digit})*
integer               {dec_digit}+
label                 [1-9]{dec_digit}*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

double                (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))

x                     {double}
y                     {double}
z                     {double}
value                 {double}

node                  ^{space}"N"{cspace}
element               ^{space}"EN"{cspace}
bface                 ^{space}"SFE"{cspace}


%%

%{
    labelList labels(8);
%}


 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

{node}{label}{cspace}{x}{cspace}{y}{cspace}{z}{space}\n {
        IStringStream nodeStream(YYText());
        char tag, c;
        label nodei;
        point node;
        nodeStream
            >> tag
            >> c >> nodei
            >> c >> node.x()
            >> c >> node.y()
            >> c >> node.z();

        if (nodei > maxNodei) maxNodei = nodei;

        slPointMap.append(nodei);
        slPoints.append(node);
    }


{element}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{space}\n {
        IStringStream elementStream(YYText());
        char tag, c;
        label celli;
        elementStream
            >> tag >> tag
            >> c >> celli
            >> c >> labels[0]
            >> c >> labels[1]
            >> c >> labels[2]
            >> c >> labels[3]
            >> c >> labels[4]
            >> c >> labels[5]
            >> c >> labels[6]
            >> c >> labels[7];

        if (celli > maxCelli) maxCelli = celli;

        slCellMap.append(celli);
        slCellLabels.append(new labelList(labels));
    }


{bface}{label}{cspace}{label}{cspace}{identifier}{cspace}{integer}{cspace}{value}{space}\n {
        IStringStream bfaceStream(YYText());
        char tag, c;
        label elementi;
        label facei;
        scalar indexValue, unknown;
        bfaceStream
            >> tag >> tag >> tag
            >> c >> elementi
            >> c >> facei
            >> c >> tag >> tag >> tag >> tag
            >> c >> unknown
            >> c >> indexValue;

        label patchi = label(indexValue);

        if (patchi > slPatchCells.size())
        {
            slPatchCells.setSize(patchi);

            forAll(slPatchCells, i)
            {
                if (!slPatchCells(i))
                {
                    slPatchCells.set(i, new SLList<label>);
                }
            }
        }

        if (patchi > slPatchCellFaces.size())
        {
            slPatchCellFaces.setSize(patchi);

            forAll(slPatchCells, i)
            {
                if (!slPatchCellFaces(i))
                {
                    slPatchCellFaces.set(i, new SLList<label>);
                }
            }
        }

        slPatchCells[patchi-1].append(elementi);
        slPatchCellFaces[patchi-1].append(facei);
    }


 /* ------------------------------------------------------------------------- *\
    ------ Ignore remaining space and \n s.  Any other characters are errors.
 \* ------------------------------------------------------------------------- */

.|\n {}


 /* ------------------------------------------------------------------------- *\
    ------ On EOF return to previous file, if none exists terminate.
 \* ------------------------------------------------------------------------- */

<<EOF>> {
            yyterminate();
    }
%%


#include "fileName.H"
#include <fstream>
using std::ifstream;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("ANSYS input file");
    argList::validOptions.insert("scale", "scale factor");

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    scalar scaleFactor = 1.0;
    args.optionReadIfPresent("scale", scaleFactor);

#   include "createTime.H"

    fileName ansysFile(args.additionalArgs()[0]);
    ifstream ansysStream(ansysFile.c_str());

    if (!ansysStream)
    {
        FatalErrorIn("ansysToFoam::main(int argc, char *argv[])")
            << args.executable()
            << ": file " << ansysFile << " not found"
            << exit(FatalError);
    }

    yyFlexLexer lexer(&ansysStream);
    while(lexer.yylex() != 0)
    {}

    Info << "Creating points" << endl;

    pointField points(slPoints.size());

    label i=0;
    for
    (
        SLList<point>::iterator pointIter = slPoints.begin();
        pointIter != slPoints.end();
        ++pointIter
    )
    {
        // Scale points for the given scale factor
        points[i++] = scaleFactor * pointIter();
    }


    labelList pointMap(maxNodei+1);

    i=0;
    for
    (
        SLList<label>::iterator pointMapIter = slPointMap.begin();
        pointMapIter != slPointMap.end();
        ++pointMapIter
    )
    {
        pointMap[pointMapIter()] = i++;
    }

    Info << "Creating cells" << endl;

    labelList cellMap(maxCelli+1);

    i=0;
    for
    (
        SLList<label>::iterator cellMapIter = slCellMap.begin();
        cellMapIter != slCellMap.end();
        ++cellMapIter
    )
    {
        cellMap[cellMapIter()] = i++;
    }


    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& tet = *(cellModeller::lookup("tet"));

    labelList labelsHex(8);
    labelList labelsPrism(6);
    labelList labelsPyramid(5);
    labelList labelsTet(4);

    cellShapeList cellShapes(slCellLabels.size());
    label nCells = 0;

    for
    (
        SLPtrList<labelList>::iterator cellIter = slCellLabels.begin();
        cellIter != slCellLabels.end();
        ++cellIter
    )
    {
        if      // Tetrahedron
        (
            cellIter()[2] == cellIter()[3]
         && cellIter()[4] == cellIter()[5]
         && cellIter()[5] == cellIter()[6]
         && cellIter()[6] == cellIter()[7]
        )
        {
            labelsTet[0] = pointMap[cellIter()[0] ];
            labelsTet[1] = pointMap[cellIter()[1] ];
            labelsTet[2] = pointMap[cellIter()[2] ];
            labelsTet[3] = pointMap[cellIter()[4] ];

            cellShapes[nCells++] = cellShape(tet, labelsTet);
        }

        else if // Square-based pyramid
        (
            cellIter()[4] == cellIter()[5]
         && cellIter()[5] == cellIter()[6]
         && cellIter()[6] == cellIter()[7]
        )
        {
            labelsPyramid[0] = pointMap[cellIter()[0] ];
            labelsPyramid[1] = pointMap[cellIter()[1] ];
            labelsPyramid[2] = pointMap[cellIter()[2] ];
            labelsPyramid[3] = pointMap[cellIter()[3] ];
            labelsPyramid[4] = pointMap[cellIter()[4] ];

            cellShapes[nCells++] = cellShape(pyr, labelsPyramid);
        }

        else if // Triangular prism
        (
            cellIter()[2] == cellIter()[3]
         && cellIter()[6] == cellIter()[7]
        )
        {
            labelsPrism[0] = pointMap[cellIter()[0] ];
            labelsPrism[1] = pointMap[cellIter()[1] ];
            labelsPrism[2] = pointMap[cellIter()[2] ];
            labelsPrism[3] = pointMap[cellIter()[4] ];
            labelsPrism[4] = pointMap[cellIter()[5] ];
            labelsPrism[5] = pointMap[cellIter()[6] ];

            cellShapes[nCells++] = cellShape(prism, labelsPrism);
        }

        else // Hex
        {
            labelsHex[0] = pointMap[cellIter()[0] ];
            labelsHex[1] = pointMap[cellIter()[1] ];
            labelsHex[2] = pointMap[cellIter()[2] ];
            labelsHex[3] = pointMap[cellIter()[3] ];
            labelsHex[4] = pointMap[cellIter()[4] ];
            labelsHex[5] = pointMap[cellIter()[5] ];
            labelsHex[6] = pointMap[cellIter()[6] ];
            labelsHex[7] = pointMap[cellIter()[7] ];

            cellShapes[nCells++] = cellShape(hex, labelsHex);
        }
    }

    // Warning: tet face order has changed between version 1.9.6 and 2.0
    //
    label faceIndex[7][6] =
    {
        {-1, -1, -1, -1, -1, -1}, // 0
        {-1, -1, -1, -1, -1, -1}, // 1
        {-1, -1, -1, -1, -1, -1}, // 2
        {-1, -1, -1, -1, -1, -1}, // 3
        {-1,  2,  0,  3,  1, -1}, // tet (version 2.0)
        { 0,  4,  3, -1,  2,  1}, // prism
        { 4,  2,  1,  3,  0,  5}, // hex
    };

    Info << "Creating boundary patches" << endl;

    faceListList boundary(slPatchCells.size());
    wordList patchNames(slPatchCells.size());

    forAll(slPatchCells, patchI)
    {
        SLList<face> patchFaces;

        SLList<label>::iterator cellIter(slPatchCells[patchI].begin());
        SLList<label>::iterator faceIter(slPatchCellFaces[patchI].begin());

        for
        (
            ;
            cellIter != slPatchCells[patchI].end()
         && faceIter != slPatchCellFaces[patchI].end();
            ++cellIter, ++faceIter
        )
        {
            patchFaces.append
            (
                cellShapes[cellMap[cellIter()] ].faces()
                [
                    faceIndex
                        [cellShapes[cellMap[cellIter()] ].nFaces()]
                        [faceIter()-1]
                ]
            );
        }

        boundary[patchI] = patchFaces;
        patchNames[patchI] = word("patch") + name(patchI + 1);
        Info << "Patch " << patchI << " named " << patchNames[patchI]
            << ": " << boundary[patchI].size() << " faces" << endl;
    }

    Info << "ansysToFoam: " << endl
        << "Ansys file format does not provide information about the type of "
        << "the patch (eg. wall, symmetry plane, cyclic etc)." << endl
        << "All the patches have been created "
        << "as type patch. Please reset after mesh conversion as necessary."
        << endl;

    wordList patchTypes(boundary.size(), polyPatch::typeName);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    wordList patchPhysicalTypes(boundary.size());

    preservePatchTypes
    (
        runTime,
        runTime.constant(),
        polyMesh::defaultRegion,
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );

    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferMove(points),
        cellShapes,
        boundary,
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << "Writing polyMesh" << endl;
    pShapeMesh.write();

    Info<< nl << "end" << endl;
    return 0;
}


 /* ------------------------------------------------------------------------- *\
    ------ End of ansysToFoam.L
 \* ------------------------------------------------------------------------- */
