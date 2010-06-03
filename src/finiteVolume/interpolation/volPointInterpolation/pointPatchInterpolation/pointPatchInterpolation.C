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

\*---------------------------------------------------------------------------*/

#include "pointPatchInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyFvPatch.H"
#include "demandDrivenData.H"
#include "coupledPointPatchFields.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointPatchInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void pointPatchInterpolation::makePatchPatchAddressing()
{
    if (debug)
    {
        Info<< "pointPatchInterpolation::makePatchPatchAddressing() : "
            << "constructing boundary addressing"
            << endl;
    }

    const fvBoundaryMesh& bm = fvMesh_.boundary();
    const pointBoundaryMesh& pbm = pointMesh::New(fvMesh_).boundary();

    // first count the total number of patch-patch points

    label nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if(!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
        {
            nPatchPatchPoints += bm[patchi].patch().boundaryPoints().size();
        }
    }


    // Go through all patches and mark up the external edge points
    Map<label> patchPatchPointSet(2*nPatchPatchPoints);

    patchPatchPoints_.setSize(nPatchPatchPoints);

    List<pointConstraint> patchPatchPointConstraints(nPatchPatchPoints);

    label pppi = 0;

    forAll(bm, patchi)
    {
        if(!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
        {
            const labelList& bp = bm[patchi].patch().boundaryPoints();
            const labelList& meshPoints = bm[patchi].patch().meshPoints();

            forAll(bp, pointi)
            {
                label ppp = meshPoints[bp[pointi]];

                Map<label>::iterator iter = patchPatchPointSet.find(ppp);

                if (iter == patchPatchPointSet.end())
                {
                    patchPatchPointSet.insert(ppp, pppi);
                    patchPatchPoints_[pppi] = ppp;

                    pbm[patchi].applyConstraint
                    (
                        bp[pointi],
                        patchPatchPointConstraints[pppi]
                    );
                    pppi++;
                }
                else
                {
                    pbm[patchi].applyConstraint
                    (
                        bp[pointi],
                        patchPatchPointConstraints[iter()]
                    );
                }
            }
        }
    }

    nPatchPatchPoints = pppi;
    patchPatchPoints_.setSize(nPatchPatchPoints);
    patchPatchPointConstraints.setSize(nPatchPatchPoints);

    patchPatchPointConstraintPoints_.setSize(nPatchPatchPoints);
    patchPatchPointConstraintTensors_.setSize(nPatchPatchPoints);

    label nConstraints = 0;

    forAll(patchPatchPointConstraints, i)
    {
        if (patchPatchPointConstraints[i].first() != 0)
        {
            patchPatchPointConstraintPoints_[nConstraints] = 
                patchPatchPoints_[i];

            patchPatchPointConstraintTensors_[nConstraints] = 
                patchPatchPointConstraints[i].constraintTransformation();

            nConstraints++;
        }
    }

    patchPatchPointConstraintPoints_.setSize(nConstraints);
    patchPatchPointConstraintTensors_.setSize(nConstraints);


    patchInterpolators_.clear();
    patchInterpolators_.setSize(bm.size());

    forAll(bm, patchi)
    {
        patchInterpolators_.set
        (
            patchi,
            new primitivePatchInterpolation(bm[patchi].patch())
        );
    }

    if (debug)
    {
        Info<< "pointPatchInterpolation::makePatchPatchAddressing() : "
            << "finished constructing boundary addressing"
            << endl;
    }
}


void pointPatchInterpolation::makePatchPatchWeights()
{
    if (debug)
    {
        Info<< "pointPatchInterpolation::makePatchPatchWeights() : "
            << "constructing boundary weighting factors"
            << endl;
    }

    patchPatchPointWeights_.clear();
    patchPatchPointWeights_.setSize(patchPatchPoints_.size());

    const labelListList& pf = fvMesh_.pointFaces();
    const volVectorField& centres = fvMesh_.C();
    const fvBoundaryMesh& bm = fvMesh_.boundary();

    pointScalarField sumWeights
    (
        IOobject
        (
            "sumWeights",
            fvMesh_.polyMesh::instance(),
            fvMesh_
        ),
        pointMesh::New(fvMesh_),
        dimensionedScalar("zero", dimless, 0)
    );

    forAll(patchPatchPoints_, pointi)
    {
        const label curPoint = patchPatchPoints_[pointi];
        const labelList& curFaces = pf[curPoint];

        patchPatchPointWeights_[pointi].setSize(curFaces.size());
        scalarList& pw = patchPatchPointWeights_[pointi];

        label nFacesAroundPoint = 0;

        const vector& pointLoc = fvMesh_.points()[curPoint];

        forAll(curFaces, facei)
        {
            if (!fvMesh_.isInternalFace(curFaces[facei]))
            {
                label patchi =
                    fvMesh_.boundaryMesh().whichPatch(curFaces[facei]);

                if (!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
                {
                    vector d =
                        pointLoc
                      - centres.boundaryField()[patchi]
                            [bm[patchi].patch().whichFace(curFaces[facei])];

                    pw[nFacesAroundPoint] = 1.0/(mag(d)+VSMALL);

                    nFacesAroundPoint++;
                }
            }
        }

        // Reset the sizes of the local weights
        pw.setSize(nFacesAroundPoint);

        // Collect the sum of weights for parallel correction
        sumWeights[curPoint] += sum(pw);
    }

    // Do parallel correction of weights

    // Update coupled boundaries
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


    // Re-scale the weights for the current point
    forAll(patchPatchPoints_, pointi)
    {
        scalarList& pw = patchPatchPointWeights_[pointi];
        scalar sumw = sumWeights[patchPatchPoints_[pointi]];

        forAll(pw, facei)
        {
            pw[facei] /= sumw;
        }
    }


    if (debug)
    {
        Info<< "pointPatchInterpolation::makePatchPatchWeights() : "
            << "finished constructing boundary weighting factors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

pointPatchInterpolation::pointPatchInterpolation(const fvMesh& vm)
:
    fvMesh_(vm)
{
    updateMesh();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

pointPatchInterpolation::~pointPatchInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointPatchInterpolation::updateMesh()
{
    makePatchPatchAddressing();
    makePatchPatchWeights();
}


bool pointPatchInterpolation::movePoints()
{
    forAll(patchInterpolators_, patchi)
    {
        patchInterpolators_[patchi].movePoints();
    }

    makePatchPatchWeights();

    return true;
}


// Specialisaion of applyCornerConstraints for scalars because
// no constraint need be applied
template<>
void pointPatchInterpolation::applyCornerConstraints<scalar>
(
    GeometricField<scalar, pointPatchField, pointMesh>& pf
) const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
