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

\*---------------------------------------------------------------------------*/

#include "pointMesh.H"
#include "globalMeshData.H"
#include "globalPointPatch.H"
#include "pointMeshMapper.H"
#include "pointFields.H"
#include "MapGeometricFields.H"
#include "MapPointField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMesh::mapFields(const mapPolyMesh& mpm)
{
    // Create a mapper
    const pointMeshMapper m(*this, mpm);

    MapGeometricFields<scalar, pointPatchField, pointMeshMapper, pointMesh>(m);
    MapGeometricFields<vector, pointPatchField, pointMeshMapper, pointMesh>(m);
    MapGeometricFields
    <
        sphericalTensor,
        pointPatchField,
        pointMeshMapper,
        pointMesh
    >(m);
    MapGeometricFields<symmTensor, pointPatchField, pointMeshMapper, pointMesh>
    (m);
    MapGeometricFields<tensor, pointPatchField, pointMeshMapper, pointMesh>(m);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMesh::pointMesh
(
    const polyMesh& pMesh,
    bool alwaysConstructGlobalPatch
)
:
    MeshObject<polyMesh, pointMesh>(pMesh),
    GeoMesh<polyMesh>(pMesh),
    boundary_(*this, pMesh.boundaryMesh())
{
    // Add the globalPointPatch if there are global points
    if
    (
        alwaysConstructGlobalPatch
     || GeoMesh<polyMesh>::mesh_.globalData().nGlobalPoints()
    )
    {
        boundary_.setSize(boundary_.size() + 1);

        boundary_.set
        (
            boundary_.size() - 1,
            new globalPointPatch
            (
                boundary_,
                boundary_.size() - 1
            )
        );
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


void Foam::pointMesh::movePoints(const pointField& newPoints)
{
    boundary_.movePoints(newPoints);
}


void Foam::pointMesh::updateMesh(const mapPolyMesh& mpm)
{
    boundary_.updateMesh();

    // Map all registered point fields
    mapFields(mpm);
}


// ************************************************************************* //
