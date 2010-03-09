/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "sampledTriSurfaceMesh.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "volFields.H"
#include "meshSearch.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledTriSurfaceMesh, 0);
    addToRunTimeSelectionTable
    (
        sampledSurface,
        sampledTriSurfaceMesh,
        word
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledTriSurfaceMesh::sampledTriSurfaceMesh
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName
)
:
    sampledSurface(name, mesh),
    surface_
    (
        IOobject
        (
            name,
            mesh.time().constant(), // instance
            "triSurface",       // local
            mesh,               // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    needsUpdate_(true),
    cellLabels_(0),
    pointMap_(0),
    faceMap_(0)
{}


Foam::sampledTriSurfaceMesh::sampledTriSurfaceMesh
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    surface_
    (
        IOobject
        (
            dict.lookup("surface"),
            mesh.time().constant(), // instance
            "triSurface",       // local
            mesh,               // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    needsUpdate_(true),
    cellLabels_(0),
    pointMap_(0),
    faceMap_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledTriSurfaceMesh::~sampledTriSurfaceMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledTriSurfaceMesh::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledTriSurfaceMesh::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();
    MeshStorage::clear();
    cellLabels_.clear();
    pointMap_.clear();
    faceMap_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledTriSurfaceMesh::update()
{
    if (!needsUpdate_)
    {
        return false;
    }


    // Find the cells the triangles of the surface are in.
    // Does approximation by looking at the face centres only
    const pointField& fc = surface_.faceCentres();

    cellLabels_.setSize(fc.size());

    meshSearch meshSearcher(mesh(), false);

    forAll(fc, triI)
    {
        cellLabels_[triI] = meshSearcher.findCell
        (
            fc[triI],
            -1,                 // seed
            true                // use tree
        );
    }


    boolList include(surface_.size(), true);

    label nNotFound = 0;

    forAll(cellLabels_, triI)
    {
        if (cellLabels_[triI] == -1)
        {
            include[triI] = false;
            nNotFound++;
        }
    }
    label nTotalNotFound = returnReduce(nNotFound, sumOp<label>());

    if (debug)
    {
        Pout<< "surface:" << surface_.size()
            << "  included:" << surface_.size()-nNotFound
            << "  total not found" << nTotalNotFound << endl;
    }


    // Add to master all triangles are outside all meshes.
    {
        boolList onAnyProc(include);
        Pstream::listCombineGather(onAnyProc, orEqOp<bool>());

        if (Pstream::master())
        {
            label nAdded = 0;
            forAll(onAnyProc, triI)
            {
                if (!onAnyProc[triI])
                {
                    include[triI] = true;
                    nAdded++;
                }
            }

            if (debug)
            {
                Pout<< "nAdded to master:" << nAdded << endl;
            }
        }
    }


    // Now subset the surface
    triSurface localSurface
    (
        surface_.subsetMesh
        (
            include,
            pointMap_,
            faceMap_
        )
    );

    // And convert into faces.
    faceList& faces = this->storedFaces();
    faces.setSize(localSurface.size());
    forAll(localSurface, i)
    {
        faces[i] = localSurface[i].triFaceFace();
    }

    this->storedPoints() = localSurface.points();

    if (debug)
    {
        print(Pout);
        Pout << endl;
    }

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField>
Foam::sampledTriSurfaceMesh::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledTriSurfaceMesh::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledTriSurfaceMesh::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledTriSurfaceMesh::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledTriSurfaceMesh::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledTriSurfaceMesh::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledTriSurfaceMesh::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledTriSurfaceMesh::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledTriSurfaceMesh::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledTriSurfaceMesh::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledTriSurfaceMesh::print(Ostream& os) const
{
    os  << "sampledTriSurfaceMesh: " << name() << " :"
        << "  surface:" << surface_.objectRegistry::name()
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
