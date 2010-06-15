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
    The function renumbers the addressing such that the band of the
    matrix is reduced. The algorithm uses a simple search through the
    neighbour list


\*---------------------------------------------------------------------------*/

#include "bandCompression.H"
#include "SLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Constructor from components
Foam::labelList Foam::bandCompression(const labelListList& cellCellAddressing)
{
    labelList newOrder(cellCellAddressing.size());

    // the business bit of the renumbering
    SLList<label> nextCell;

    labelList visited(cellCellAddressing.size());

    label currentCell;
    label cellInOrder = 0;

    // reset the visited cells list
    forAll (visited, cellI)
    {
        visited[cellI] = 0;
    }

    // loop over the cells
    forAll (visited, cellI)
    {
        // find the first cell that has not been visited yet
        if (visited[cellI] == 0)
        {
            currentCell = cellI;

            // use this cell as a start
            nextCell.append(currentCell);

            // loop through the nextCell list. Add the first cell into the
            // cell order if it has not already been visited and ask for its
            // neighbours. If the neighbour in question has not been visited,
            // add it to the end of the nextCell list

            while (nextCell.size())
            {
                currentCell = nextCell.removeHead();

                if (visited[currentCell] == 0)
                {
                    visited[currentCell] = 1;

                    // add into cellOrder
                    newOrder[cellInOrder] = currentCell;
                    cellInOrder++;

                    // find if the neighbours have been visited
                    const labelList& neighbours =
                        cellCellAddressing[currentCell];

                    forAll (neighbours, nI)
                    {
                        if (visited[neighbours[nI]] == 0)
                        {
                            // not visited, add to the list
                            nextCell.append(neighbours[nI]);
                        }
                    }
                }
            }
        }
    }

    return newOrder;
}


// ************************************************************************* //
