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

#include "linearUpwind.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::linearUpwind<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "linearUpwindCorrection(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr();

    const surfaceScalarField& faceFlux = this->faceFlux_;

    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    GeometricField
        <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        gradVf = gradScheme_().grad(vf);

    forAll(faceFlux, facei)
    {
        if (faceFlux[facei] > 0)
        {
            label own = owner[facei];
            sfCorr[facei] = (Cf[facei] - C[own]) & gradVf[own];
        }
        else
        {
            label nei = neighbour[facei];
            sfCorr[facei] = (Cf[facei] - C[nei]) & gradVf[nei];
        }
    }


    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        GeometricBoundaryField& bSfCorr = sfCorr.boundaryField();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const unallocLabelList& pOwner = 
                mesh.boundary()[patchi].faceCells();

            const vectorField& pCf = Cf.boundaryField()[patchi];

            const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

            Field<typename outerProduct<vector, Type>::type> pGradVfNei =
                gradVf.boundaryField()[patchi].patchNeighbourField();

            // Build the d-vectors
            vectorField pd = 
                mesh.Sf().boundaryField()[patchi]
               /(
                   mesh.magSf().boundaryField()[patchi]
                  *mesh.deltaCoeffs().boundaryField()[patchi]
                );

            if (!mesh.orthogonal())
            {
                pd -= mesh.correctionVectors().boundaryField()[patchi]
                    /mesh.deltaCoeffs().boundaryField()[patchi];
            }

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                if (pFaceFlux[facei] > 0)
                {
                    pSfCorr[facei] = (pCf[facei] - C[own]) & gradVf[own];
                }
                else
                {
                    pSfCorr[facei] = 
                        (pCf[facei] - pd[facei] - C[own]) & pGradVfNei[facei];
                }
            }
        }
    }

    return tsfCorr;
}


namespace Foam
{
    //makelimitedSurfaceInterpolationScheme(linearUpwind)
    makelimitedSurfaceInterpolationTypeScheme(linearUpwind, scalar)
    makelimitedSurfaceInterpolationTypeScheme(linearUpwind, vector)
}

// ************************************************************************* //
