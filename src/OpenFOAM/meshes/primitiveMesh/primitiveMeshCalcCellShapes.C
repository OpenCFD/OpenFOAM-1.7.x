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

#include "primitiveMesh.H"
#include "degenerateMatcher.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellShapes() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcCellShapes() : calculating cellShapes"
            << endl;
    }

    // It is an error to attempt to recalculate faceCells
    // if the pointer is already set
    if (cellShapesPtr_)
    {
        FatalErrorIn("primitiveMesh::calcCellShapes() const")
            << "cellShapes already calculated"
            << abort(FatalError);
    }
    else
    {
        cellShapesPtr_ = new cellShapeList(nCells());
        cellShapeList& cellShapes = *cellShapesPtr_;

        forAll (cellShapes, celli)
        {
            cellShapes[celli] = degenerateMatcher::match(*this, celli);
        }
    }
}


// ************************************************************************* //
