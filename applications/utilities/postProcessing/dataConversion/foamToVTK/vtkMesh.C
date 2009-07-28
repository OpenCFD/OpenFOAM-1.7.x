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

#include "vtkMesh.H"
#include "fvMeshSubset.H"
#include "Time.H"
#include "cellSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::vtkMesh::vtkMesh
(
    fvMesh& baseMesh,
    const word& setName
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    setName_(setName)
{
    if (setName.size())
    {
        // Read cellSet using whole mesh
        cellSet currentSet(baseMesh_, setName_);

        // Set current subset
        subsetter_.setLargeCellSubset(currentSet);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::polyMesh::readUpdateState Foam::vtkMesh::readUpdate()
{
    polyMesh::readUpdateState meshState = baseMesh_.readUpdate();

    if (meshState != polyMesh::UNCHANGED)
    {
        // Note: since fvMeshSubset has no movePoints() functionality reconstruct
        // the subset even if only movement.

        topoPtr_.clear();

        if (setName_.size())
        {
            Info<< "Subsetting mesh based on cellSet " << setName_ << endl;

            // Read cellSet using whole mesh
            cellSet currentSet(baseMesh_, setName_);

            subsetter_.setLargeCellSubset(currentSet);
        }
    }

    return meshState;
}


// ************************************************************************* //
