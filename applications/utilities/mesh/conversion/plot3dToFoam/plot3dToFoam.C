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
    Plot3d mesh (ascii/formatted format) converter.

    Work in progress! Handles ascii multiblock (and optionally singleBlock)
    format.
    By default expects blanking. Use -noBlank if none.
    Use -2D @a thickness if 2D.

    Niklas Nordin has experienced a problem with lefthandedness of the blocks.
    The code should detect this automatically - see hexBlock::readPoints but
    if this goes wrong just set the blockHandedness_ variable to 'right'
    always.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IFstream.H"
#include "hexBlock.H"
#include "polyMesh.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "preservePatchTypes.H"
#include "cellShape.H"
#include "cellModeller.H"
#include "mergePoints.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("PLOT3D geom file");
    argList::validOptions.insert("scale", "scale factor");
    argList::validOptions.insert("noBlank", "");
    argList::validOptions.insert("singleBlock", "");
    argList::validOptions.insert("2D", "thickness");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    scalar scaleFactor = 1.0;
    args.optionReadIfPresent("scale", scaleFactor);

    bool readBlank = !args.optionFound("noBlank");
    bool singleBlock = args.optionFound("singleBlock");
    scalar twoDThickness = -1;
    if (args.optionReadIfPresent("2D", twoDThickness))
    {
        Info<< "Reading 2D case by extruding points by " << twoDThickness
            << " in z direction." << nl << endl;
    }


#   include "createTime.H"

    IFstream plot3dFile(args.additionalArgs()[0]);

    // Read the plot3d information using a fixed format reader.
    // Comments in the file are in C++ style, so the stream parser will remove
    // them with no intervention
    label nblock;

    if (singleBlock)
    {
        nblock = 1;
    }
    else
    {
        plot3dFile >> nblock;
    }

    Info<< "Reading " << nblock << " blocks" << endl;

    PtrList<hexBlock> blocks(nblock);

    {
        label nx, ny, nz;

        forAll (blocks, blockI)
        {
            if (twoDThickness > 0)
            {
                // Fake second set of points (done in readPoints below)
                plot3dFile >> nx >> ny;
                nz = 2;
            }
            else
            {
                plot3dFile >> nx >> ny >> nz;
            }

            Info<< "block " << blockI << " nx:" << nx
                << " ny:" << ny << " nz:" << nz << endl;

            blocks.set(blockI, new hexBlock(nx, ny, nz));
        }
    }

    Info<< "Reading block points" << endl;
    label sumPoints(0);
    label nMeshCells(0);

    forAll (blocks, blockI)
    {
        Info<< "block " << blockI << ":" << nl;
        blocks[blockI].readPoints(readBlank, twoDThickness, plot3dFile);
        sumPoints += blocks[blockI].nBlockPoints();
        nMeshCells += blocks[blockI].nBlockCells();
        Info<< nl;
    }

    pointField points(sumPoints);
    labelList blockOffsets(blocks.size());
    sumPoints = 0;
    forAll (blocks, blockI)
    {
        const pointField& blockPoints = blocks[blockI].points();
        blockOffsets[blockI] = sumPoints;
        forAll (blockPoints, i)
        {
            points[sumPoints++] = blockPoints[i];
        }
    }

    // From old to new master point
    labelList oldToNew;
    pointField newPoints;

    // Merge points
    mergePoints
    (
        points,
        SMALL,
        false,
        oldToNew,
        newPoints
    );

    Info<< "Merged points within " << SMALL << " distance. Merged from "
        << oldToNew.size() << " down to " << newPoints.size()
        << " points." << endl;

    // Scale the points
    if (scaleFactor > 1.0 + SMALL || scaleFactor < 1.0 - SMALL)
    {
        newPoints *= scaleFactor;
    }

    Info<< "Creating cells" << endl;

    cellShapeList cellShapes(nMeshCells);

    const cellModel& hex = *(cellModeller::lookup("hex"));

    label nCreatedCells = 0;

    forAll (blocks, blockI)
    {
        labelListList curBlockCells = blocks[blockI].blockCells();

        forAll (curBlockCells, blockCellI)
        {
            labelList cellPoints(curBlockCells[blockCellI].size());

            forAll (cellPoints, pointI)
            {
                cellPoints[pointI] =
                    oldToNew
                    [
                        curBlockCells[blockCellI][pointI]
                      + blockOffsets[blockI]
                    ];
            }

            // Do automatic collapse from hex.
            cellShapes[nCreatedCells] = cellShape(hex, cellPoints, true);

            nCreatedCells++;
        }
    }

    Info<< "Creating boundary patches" << endl;

    faceListList boundary(0);
    wordList patchNames(0);
    wordList patchTypes(0);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = wallPolyPatch::typeName;
    wordList patchPhysicalTypes(0);

    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferMove(newPoints),
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

    Info<< "Writing polyMesh" << endl;
    pShapeMesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
