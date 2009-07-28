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

#include "refinementParameters.H"
#include "mathematicalConstants.H"
#include "polyMesh.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::refinementParameters::refinementParameters
(
    const dictionary& dict,
    const label dummy
)
:
    maxGlobalCells_(readLabel(dict.lookup("cellLimit"))),
    maxLocalCells_(readLabel(dict.lookup("procCellLimit"))),
    minRefineCells_(readLabel(dict.lookup("minimumRefine"))),
    curvature_(readScalar(dict.lookup("curvature"))),
    nBufferLayers_(readLabel(dict.lookup("nBufferLayers"))),
    keepPoints_(dict.lookup("keepPoints"))
{}


Foam::refinementParameters::refinementParameters(const dictionary& dict)
:
    maxGlobalCells_(readLabel(dict.lookup("maxGlobalCells"))),
    maxLocalCells_(readLabel(dict.lookup("maxLocalCells"))),
    minRefineCells_(readLabel(dict.lookup("minRefinementCells"))),
    nBufferLayers_(readLabel(dict.lookup("nCellsBetweenLevels"))),
    keepPoints_(pointField(1, dict.lookup("locationInMesh")))
{
    scalar featAngle(readScalar(dict.lookup("resolveFeatureAngle")));

    if (featAngle < 0 || featAngle > 180)
    {
        curvature_ = -GREAT;
    }
    else
    {
        curvature_ = Foam::cos(featAngle*mathematicalConstant::pi/180.0);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::refinementParameters::findCells(const polyMesh& mesh)
 const
{
    // Global calculation engine
    globalIndex globalCells(mesh.nCells());

    // Cell label per point
    labelList cellLabels(keepPoints_.size());

    forAll(keepPoints_, i)
    {
        const point& keepPoint = keepPoints_[i];

        label localCellI = mesh.findCell(keepPoint);

        label globalCellI = -1;

        if (localCellI != -1)
        {
            Pout<< "Found point " << keepPoint << " in cell " << localCellI
                << " on processor " << Pstream::myProcNo() << endl;
            globalCellI = globalCells.toGlobal(localCellI);
        }

        reduce(globalCellI, maxOp<label>());

        if (globalCellI == -1)
        {
            FatalErrorIn
            (
                "refinementParameters::findCells(const polyMesh&) const"
            )   << "Point " << keepPoint
                << " is not inside the mesh or on a face or edge." << nl
                << "Bounding box of the mesh:" << mesh.bounds()
                << exit(FatalError);
        }

        if (globalCells.isLocal(globalCellI))
        {
            cellLabels[i] = localCellI;
        }
        else
        {
            cellLabels[i] = -1;
        }
    }
    return cellLabels;
}


// ************************************************************************* //
