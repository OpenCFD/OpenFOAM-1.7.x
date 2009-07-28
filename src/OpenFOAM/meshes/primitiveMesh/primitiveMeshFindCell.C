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

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "cell.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Is the point in the cell bounding box
bool Foam::primitiveMesh::pointInCellBB(const point& p, label celli) const
{
    const pointField& points = this->points();
    const faceList& f = faces();
    const vectorField& centres = cellCentres();
    const cellList& cf = cells();

    labelList cellVertices = cf[celli].labels(f);

    vector bbmax = -GREAT*vector::one;
    vector bbmin = GREAT*vector::one;

    forAll (cellVertices, vertexI)
    {
        bbmax = max(bbmax, points[cellVertices[vertexI]]);
        bbmin = min(bbmin, points[cellVertices[vertexI]]);
    }

    scalar distance = mag(centres[celli] - p);

    if ((distance - mag(bbmax - bbmin)) < SMALL)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Is the point in the cell
bool Foam::primitiveMesh::pointInCell(const point& p, label celli) const
{
    const labelList& f = cells()[celli];
    const labelList& owner = this->faceOwner();
    const vectorField& cf = faceCentres();
    const vectorField& Sf = faceAreas();

    bool inCell = true;

    forAll(f, facei)
    {
        label nFace = f[facei];
        vector proj = p - cf[nFace];
        vector normal = Sf[nFace];
        if (owner[nFace] != celli)
        {
            normal = -normal;
        }
        inCell = inCell && ((normal & proj) <= 0);
    }

    return inCell;
}


// Find the cell with the nearest cell centre
Foam::label Foam::primitiveMesh::findNearestCell(const point& location) const
{
    const vectorField& centres = cellCentres();

    label nearestCelli = 0;
    scalar minProximity = magSqr(centres[0] - location);

    for (label celli = 1; celli < centres.size(); celli++)
    {
        scalar proximity = magSqr(centres[celli] - location);

        if (proximity < minProximity)
        {
            nearestCelli = celli;
            minProximity = proximity;
        }
    }

    return nearestCelli;
}


// Find cell enclosing this location
Foam::label Foam::primitiveMesh::findCell(const point& location) const
{
    if (nCells() == 0)
    {
        return -1;
    }

    // Find the nearest cell centre to this location
    label celli = findNearestCell(location);

    // If point is in the nearest cell return
    if (pointInCell(location, celli))
    {
        return celli;
    }
    else // point is not in the nearest cell so search all cells
    {
        bool cellFound = false;
        label n = 0;

        while ((!cellFound) && (n < nCells()))
        {
            if (pointInCell(location, n))
            {
                cellFound = true;
                celli = n;
            }
            else
            {
                n++;
            }
        }
        if (cellFound)
        {
            return celli;
        }
        else
        {
            return -1;
        }
    }
}


// ************************************************************************* //
