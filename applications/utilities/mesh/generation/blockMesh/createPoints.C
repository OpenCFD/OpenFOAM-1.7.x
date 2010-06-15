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


\*---------------------------------------------------------------------------*/

#include "error.H"
#include "blockMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointField Foam::blockMesh::createPoints(const dictionary& dict)
{
    blockMesh& blocks = *this;

    scalar scaleFactor = 1.0;

    // optional 'convertToMeters' (or 'scale'?)
    if (!dict.readIfPresent("convertToMeters", scaleFactor))
    {
        dict.readIfPresent("scale", scaleFactor);
    }

    Info<< nl << "Creating points with scale " << scaleFactor << endl;

    pointField points(nPoints_);

    forAll(blocks, blockLabel)
    {
        const pointField& blockPoints = blocks[blockLabel].points();

        forAll(blockPoints, blockPointLabel)
        {
            points
            [
                mergeList_
                [
                    blockPointLabel
                  + blockOffsets_[blockLabel]
                ]
            ] = scaleFactor * blockPoints[blockPointLabel];
        }
    }

    return points;
}

// ************************************************************************* //
