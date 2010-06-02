/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Create intermediate mesh files from PROSTAR files

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void starMesh::writeMesh()
{
    if (isShapeMesh_)
    {
        Info << "This is a shapeMesh." << endl;

        Info << "Default patch type set to empty" << endl;

        clearExtraStorage();

        polyMesh pShapeMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime_.constant(),
                runTime_
            ),
            xferCopy(points_),           // we could probably re-use the data
            cellShapes_,
            boundary_,
            patchNames_,
            patchTypes_,
            defaultFacesName_,
            defaultFacesType_,
            patchPhysicalTypes_
        );

        Info << "Writing polyMesh" << endl;
        pShapeMesh.write();
    }
    else
    {
        // This is a polyMesh.

        createPolyMeshData();

        Info << "This is a polyMesh" << endl;

        clearExtraStorage();

        polyMesh pMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime_.constant(),
                runTime_
            ),
            xferCopy(points_),           // we could probably re-use the data
            xferCopy(meshFaces_),
            xferCopy(cellPolys_)
        );

        // adding patches also checks the mesh
        pMesh.addPatches(polyBoundaryPatches(pMesh));

        Info << "Writing polyMesh" << endl;
        pMesh.write();
    }
}


// ************************************************************************* //
