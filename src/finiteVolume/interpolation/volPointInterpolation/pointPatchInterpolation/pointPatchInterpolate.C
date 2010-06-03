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
#include "volFields.H"
#include "pointFields.H"
#include "emptyFvPatch.H"
#include "valuePointPatchField.H"
#include "coupledPointPatchField.H"
#include "coupledFacePointPatch.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void pointPatchInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    bool overrideFixedValue
) const
{
    if (debug)
    {
        Info<< "pointPatchInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    // Interpolate patch values: over-ride the internal values for the points
    // on the patch with the interpolated point values from the faces of the
    // patch

    const fvBoundaryMesh& bm = fvMesh_.boundary();
    const pointBoundaryMesh& pbm = pointMesh::New(fvMesh_).boundary();

    forAll(bm, patchi)
    {
        if (!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
        {
            pointPatchField<Type>& ppf = pf.boundaryField()[patchi];

            // Only map the values corresponding to the points associated with
            // faces, not "lone" points due to decomposition
            ppf.setInInternalField
            (
                pf.internalField(),
                patchInterpolators_[patchi]
               .faceToPointInterpolate(vf.boundaryField()[patchi])(),
                bm[patchi].patch().meshPoints()
            );

            if
            (
                overrideFixedValue
             && isA<valuePointPatchField<Type> >(ppf)
            )
            {
                refCast<valuePointPatchField<Type> >(ppf) = ppf;
            }
        }
        else if (bm[patchi].coupled())
        {
            // Initialise the "lone" points on the coupled patch to zero,
            // these values are obtained from the couple-transfer

            const labelList& loneMeshPoints =
                refCast<const coupledFacePointPatch>(pbm[patchi])
               .loneMeshPoints();

            forAll(loneMeshPoints, i)
            {
                pf[loneMeshPoints[i]] = pTraits<Type>::zero;
            }
        }

    }


    // Correct patch-patch boundary points by interpolation "around" corners
    const labelListList& PointFaces = fvMesh_.pointFaces();

    forAll(patchPatchPoints_, pointi)
    {
        const label curPoint = patchPatchPoints_[pointi];
        const labelList& curFaces = PointFaces[curPoint];

        label fI = 0;

        // Reset the boundary value before accumulation
        pf[curPoint] = pTraits<Type>::zero;

        // Go through all the faces
        forAll(curFaces, facei)
        {
            if (!fvMesh_.isInternalFace(curFaces[facei]))
            {
                label patchi =
                    fvMesh_.boundaryMesh().whichPatch(curFaces[facei]);

                if (!isA<emptyFvPatch>(bm[patchi]) && !bm[patchi].coupled())
                {
                    label faceInPatchi =
                        bm[patchi].patch().whichFace(curFaces[facei]);

                    pf[curPoint] +=
                        patchPatchPointWeights_[pointi][fI]
                       *vf.boundaryField()[patchi][faceInPatchi];

                    fI++;
                }
            }
        }
    }

    // Update coupled boundaries
    forAll(pf.boundaryField(), patchi)
    {
        if (pf.boundaryField()[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type> >(pf.boundaryField()[patchi])
                .initSwapAdd(pf.internalField());
        }
    }

    forAll(pf.boundaryField(), patchi)
    {
        if (pf.boundaryField()[patchi].coupled())
        {
            refCast<coupledPointPatchField<Type> >(pf.boundaryField()[patchi])
                .swapAdd(pf.internalField());
        }
    }


    // Override constrained pointPatchField types with the constraint value.
    // This relys on only constrained pointPatchField implementing the evaluate
    // function
    pf.correctBoundaryConditions();


    // Apply multiple constraints on edge/corner points
    applyCornerConstraints(pf);


    if (debug)
    {
        Info<< "pointPatchInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "finished interpolating field from cells to points"
            << endl;
    }
}


template<class Type>
void pointPatchInterpolation::applyCornerConstraints
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    forAll(patchPatchPointConstraintPoints_, pointi)
    {
        pf[patchPatchPointConstraintPoints_[pointi]] = transform
        (
            patchPatchPointConstraintTensors_[pointi],
            pf[patchPatchPointConstraintPoints_[pointi]]
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
