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

#include "volPointInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "demandDrivenData.H"
#include "coupledPointPatchFields.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void volPointInterpolation::makeWeights()
{
    if (debug)
    {
        Info<< "volPointInterpolation::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    const pointField& points = mesh().points();
    const labelListList& pointCells = mesh().pointCells();
    const vectorField& cellCentres = mesh().cellCentres();

    // Allocate storage for weighting factors
    pointWeights_.clear();
    pointWeights_.setSize(points.size());

    forAll(pointWeights_, pointi)
    {
        pointWeights_[pointi].setSize(pointCells[pointi].size());
    }

    pointScalarField sumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh().polyMesh::instance(),
            mesh()
        ),
        pointMesh::New(mesh()),
        dimensionedScalar("zero", dimless, 0)
    );

    // Calculate inverse distances between cell centres and points
    // and store in weighting factor array
    forAll(points, pointi)
    {
        scalarList& pw = pointWeights_[pointi];
        const labelList& pcp = pointCells[pointi];

        forAll(pcp, pointCelli)
        {
            pw[pointCelli] = 
                1.0/mag(points[pointi] - cellCentres[pcp[pointCelli]]);

            sumWeights[pointi] += pw[pointCelli];
        }
    }

    forAll(sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            refCast<coupledPointPatchScalarField>
                (sumWeights.boundaryField()[patchi]).initSwapAdd
                (
                    sumWeights.internalField()
                );
        }
    }

    forAll(sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            refCast<coupledPointPatchScalarField>
            (sumWeights.boundaryField()[patchi]).swapAdd
            (
                sumWeights.internalField()
            );
        }
    }

    forAll(points, pointi)
    {
        scalarList& pw = pointWeights_[pointi];

        forAll(pw, pointCelli)
        {
            pw[pointCelli] /= sumWeights[pointi];
        }
    }

    if (debug)
    {
        Info<< "volPointInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

volPointInterpolation::volPointInterpolation(const fvMesh& vm)
:
    MeshObject<fvMesh, volPointInterpolation>(vm),
    boundaryInterpolator_(vm)
{
    updateMesh();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

volPointInterpolation::~volPointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void volPointInterpolation::updateMesh()
{
    makeWeights();
    boundaryInterpolator_.updateMesh();
}


bool volPointInterpolation::movePoints()
{
    makeWeights();
    boundaryInterpolator_.movePoints();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
