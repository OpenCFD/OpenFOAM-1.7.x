/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "sampledTriSurfaceMesh.H"
#include "treeDataPoint.H"
#include "meshSearch.H"
#include "Tuple2.H"
#include "globalIndex.H"

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

    //- Private class for finding nearest
    //  - global index
    //  - sqr(distance)
    typedef Tuple2<scalar, label> nearInfo;

    class nearestEqOp
    {

    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first() < x.first())
            {
                x = y;
            }
        }
    };
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
            "triSurface",           // local
            mesh,                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    needsUpdate_(true),
    cellLabels_(0),
    pointToFace_(0)
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
            "triSurface",           // local
            mesh,                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    needsUpdate_(true),
    cellLabels_(0),
    pointToFace_(0)
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
    pointToFace_.clear();

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

    meshSearch meshSearcher(mesh(), false);

    const indexedOctree<treeDataPoint>& cellCentreTree =
        meshSearcher.cellCentreTree();


    // Global numbering for cells - only used to uniquely identify local cells.
    globalIndex globalCells(mesh().nCells());
    List<nearInfo> nearest(fc.size());
    forAll(nearest, i)
    {
        nearest[i].first() = GREAT;
        nearest[i].second() = labelMax;
    }

    // Search triangles using nearest on local mesh
    forAll(fc, triI)
    {
        pointIndexHit nearInfo = cellCentreTree.findNearest
        (
            fc[triI],
            sqr(GREAT)
        );
        if (nearInfo.hit())
        {
            nearest[triI].first() = magSqr(nearInfo.hitPoint()-fc[triI]);
            nearest[triI].second() = globalCells.toGlobal(nearInfo.index());
        }
    }

    // See which processor has the nearest.
    Pstream::listCombineGather(nearest, nearestEqOp());
    Pstream::listCombineScatter(nearest);

    boolList include(surface_.size(), false);

    cellLabels_.setSize(fc.size());
    cellLabels_ = -1;

    label nFound = 0;
    forAll(nearest, triI)
    {
        if (nearest[triI].second() == labelMax)
        {
            // Not found on any processor. How to map?
        }
        else if (globalCells.isLocal(nearest[triI].second()))
        {
            cellLabels_[triI] = globalCells.toLocal(nearest[triI].second());

            include[triI] = true;
            nFound++;
        }
    }


    if (debug)
    {
        Pout<< "Local out of faces:" << cellLabels_.size()
            << " keeping:" << nFound << endl;
    }

    // Now subset the surface. Do not use triSurface::subsetMesh since requires
    // original surface to be in compact numbering.

    const triSurface& s = surface_;

    // Compact to original triangle
    labelList faceMap(s.size());
    // Compact to original points
    labelList pointMap(s.points().size());
    // From original point to compact points
    labelList reversePointMap(s.points().size(), -1);

    {
        label newPointI = 0;
        label newTriI = 0;

        forAll(s, triI)
        {
            if (include[triI])
            {
                faceMap[newTriI++] = triI;

                const labelledTri& f = s[triI];
                forAll(f, fp)
                {
                    if (reversePointMap[f[fp]] == -1)
                    {
                        pointMap[newPointI] = f[fp];
                        reversePointMap[f[fp]] = newPointI++;
                    }
                }
            }
        }
        faceMap.setSize(newTriI);
        pointMap.setSize(newPointI);
    }

    // Subset cellLabels
    cellLabels_ = UIndirectList<label>(cellLabels_, faceMap)();

    // Store any face per point
    pointToFace_.setSize(pointMap.size());

    // Create faces and points for subsetted surface
    faceList& faces = this->storedFaces();
    faces.setSize(faceMap.size());
    forAll(faceMap, i)
    {
        const triFace& f = s[faceMap[i]];
        triFace newF
        (
            reversePointMap[f[0]],
            reversePointMap[f[1]],
            reversePointMap[f[2]]
        );
        faces[i] = newF.triFaceFace();

        forAll(newF, fp)
        {
            pointToFace_[newF[fp]] = i;
        }
    }

    this->storedPoints() = pointField(s.points(), pointMap);

    if (debug)
    {
        print(Pout);
        Pout<< endl;
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
