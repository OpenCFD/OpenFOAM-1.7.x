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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IStringStream.H"

#include "myBoundBox.H"
#include "myBoundBoxList.H"
#include "octree.H"
#include "octreeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:



int main(int argc, char *argv[])
{
    argList::validOptions.insert("x1", "X1");
    argList::validOptions.insert("y1", "Y1");
    argList::validOptions.insert("z1", "Z1");

    argList::validOptions.insert("x2", "X2");
    argList::validOptions.insert("y2", "Y2");
    argList::validOptions.insert("z2", "Z2");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Calculate BB of all cells


    myBoundBoxList allBb(mesh.nCells());

    const pointField& allPoints = mesh.points();

     vectorField bbMin(mesh.nCells());
     bbMin = vector(GREAT, GREAT, GREAT);
     vectorField bbMax(mesh.nCells());
     bbMax = vector(-GREAT, -GREAT, -GREAT);

     const labelListList& pCells = mesh.pointCells();

     forAll(pCells, pointi)
     {
         const point& vertCoord = allPoints[pointi];

         const labelList& cells = pCells[pointi];

         forAll(cells, celli)
         {
             label cellNum = cells[celli];

             bbMin[cellNum].x() = min(bbMin[cellNum].x(), vertCoord.x());
             bbMin[cellNum].y() = min(bbMin[cellNum].y(), vertCoord.y());
             bbMin[cellNum].z() = min(bbMin[cellNum].z(), vertCoord.z());

             bbMax[cellNum].x() = max(bbMax[cellNum].x(), vertCoord.x());
             bbMax[cellNum].y() = max(bbMax[cellNum].y(), vertCoord.y());
             bbMax[cellNum].z() = max(bbMax[cellNum].z(), vertCoord.z());
         }
     }


    forAll(allBb, celli)
    {
        allBb[celli] = myBoundBox(bbMin[celli], bbMax[celli]);
    }


    myBoundBox meshBb(allPoints);

    scalar typDim = meshBb.minDim()/111;

    myBoundBox shiftedBb
    (
        meshBb.min(),
        point
        (
            meshBb.max().x() + typDim,
            meshBb.max().y() + typDim,
            meshBb.max().z() + typDim
        )
    );


    Info<< "Mesh" << endl;
    Info<< "   bounding box     :" << shiftedBb << endl;
    Info<< "   typical dimension:" << shiftedBb.typDim() << endl;


    /*
     * Now we have allBb and shiftedBb
     */



    // Construct table of subset of cells

    labelList cellIndices(10);

    cellIndices[0] = 1433;
    cellIndices[1] = 1434;
    cellIndices[2] = 1435;
    cellIndices[3] = 1436;
    cellIndices[4] = 1437;
    cellIndices[5] = 1438;
    cellIndices[6] = 1439;
    cellIndices[7] = 1440;
    cellIndices[8] = 1441;
    cellIndices[9] = 1442;

    // Get the corresponding bounding boxes

    forAll(cellIndices, i)
    {
        allBb[i] = allBb[cellIndices[i]];
    }
    allBb.setSize(cellIndices.size());
    


    // Wrap indices and mesh information into helper object
    octreeData shapes(mesh, cellIndices);

    octree oc
    (
        shiftedBb,  // overall bounding box
        shapes,     // all information needed to do checks on cells
        allBb,      // bounding boxes of cells
        10.0        // maximum ratio of cubes v.s. cells
    );

//    scalar x1(readScalar(IStringStream(args.options()["x1"])()));
//    scalar y1(readScalar(IStringStream(args.options()["y1"])()));
//    scalar z1(readScalar(IStringStream(args.options()["z1"])()));

//    scalar x2(readScalar(IStringStream(args.options()["x2"])()));
//    scalar y2(readScalar(IStringStream(args.options()["y2"])()));
//    scalar z2(readScalar(IStringStream(args.options()["z2"])()));



    label nFound = 0;

    scalar x = -5.0;
    for(int i = 0; i < 100; i++)
    {
        scalar y = -7.0;
        for(int j = 0; j < 10; j++)
        {
            scalar z = -12.0;
            for (int k = 0; k < 10; k++)
            {
                point sample(x, y, z);

                label index = oc.find(sample);

                // Convert index into shapes back into cellindex.
                label cell;
                if (index != -1)
                {
                    cell = cellIndices[index];
                }
                else
                {
                    cell = -1;
                }
                Info<< "Point:" << sample
                    << " is in cell " << cell << "(octree)  "
                    << mesh.findCell(sample) << "(linear)"
                    << endl;

                z += 1.2;
            }
            y += 0.9;
        }
        x += 0.1;
    }


    Info<< "nFound=" << nFound << endl;

    Info<< "End\n" << endl;


    Info<< "Statistics:" << endl
        << "  nCells   :" << allBb.size() << endl
        << "  nNodes   :" << oc.nNodes() << endl
        << "  nLeaves  :" << oc.nLeaves() << endl
        << "  nEntries :" << oc.nEntries() << endl
        << "  Cells per leaf :"
        << oc.nEntries()/oc.nLeaves()
        << endl
        << "  Every cell in  :"
        << oc.nEntries()/allBb.size() << " cubes"
        << endl;

    return 0;
}


// ************************************************************************* //
