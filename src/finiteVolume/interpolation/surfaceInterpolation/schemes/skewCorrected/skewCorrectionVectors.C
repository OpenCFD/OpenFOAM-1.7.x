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

#include "skewCorrectionVectors.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(skewCorrectionVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::skewCorrectionVectors::skewCorrectionVectors(const fvMesh& mesh)
:
    MeshObject<fvMesh, skewCorrectionVectors>(mesh),
    skew_(true),
    skewCorrectionVectors_(NULL)
{}


Foam::skewCorrectionVectors::~skewCorrectionVectors()
{
    deleteDemandDrivenData(skewCorrectionVectors_);
}


void Foam::skewCorrectionVectors::makeSkewCorrectionVectors() const
{
    if (debug)
    {
        Info<< "surfaceInterpolation::makeSkewCorrectionVectors() : "
            << "Constructing skew correction vectors"
            << endl;
    }

    skewCorrectionVectors_ = new surfaceVectorField
    (
        IOobject
        (
            "skewCorrectionVectors",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimless
    );
    surfaceVectorField& SkewCorrVecs = *skewCorrectionVectors_;

    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const unallocLabelList& owner = mesh_.owner();

    // Build the d-vectors
    surfaceVectorField d = Sf/(mesh_.magSf()*mesh_.deltaCoeffs());

    if (!mesh_.orthogonal())
    {
        d -= mesh_.correctionVectors()/mesh_.deltaCoeffs();
    }

    forAll(owner, faceI)
    {
        vector Cpf = Cf[faceI] - C[owner[faceI]];

        SkewCorrVecs[faceI] =
            Cpf - ((Sf[faceI] & Cpf)/(Sf[faceI] & d[faceI]))*d[faceI];
    }


    forAll(SkewCorrVecs.boundaryField(), patchI)
    {
        fvsPatchVectorField& patchSkewCorrVecs =
            SkewCorrVecs.boundaryField()[patchI];

        if (!patchSkewCorrVecs.coupled())
        {
            patchSkewCorrVecs = vector::zero;
        }
        else
        {
            const fvPatch& p = patchSkewCorrVecs.patch();
            const unallocLabelList& faceCells = p.faceCells();
            const vectorField& patchFaceCentres = Cf.boundaryField()[patchI];
            const vectorField& patchSf = Sf.boundaryField()[patchI];
            const vectorField& patchD = d.boundaryField()[patchI];

            forAll(p, patchFaceI)
            {
                vector Cpf =
                    patchFaceCentres[patchFaceI] - C[faceCells[patchFaceI]];

                patchSkewCorrVecs[patchFaceI] =
                    Cpf
                  - (
                        (patchSf[patchFaceI] & Cpf)/
                        (patchSf[patchFaceI] & patchD[patchFaceI])
                    )*patchD[patchFaceI];
            }
        }
    }

    scalar skewCoeff = 0.0;

    if (Sf.internalField().size())
    {
        skewCoeff = max(mag(SkewCorrVecs)/mag(d)).value();
    }

    if (debug)
    {
        Info<< "surfaceInterpolation::makeSkewCorrectionVectors() : "
            << "skew coefficient = " << skewCoeff << endl;
    }

    //skewCoeff = 0.0;

    if (skewCoeff < 1e-5)
    {
        skew_ = false;
        deleteDemandDrivenData(skewCorrectionVectors_);
    }
    else
    {
        skew_ = true;
    }

    if (debug)
    {
        Info<< "surfaceInterpolation::makeSkewCorrectionVectors() : "
            << "Finished constructing skew correction vectors"
            << endl;
    }
}


bool Foam::skewCorrectionVectors::skew() const
{
    if (skew_ == true && !skewCorrectionVectors_)
    {
        makeSkewCorrectionVectors();
    }

    return skew_;
}


const Foam::surfaceVectorField& Foam::skewCorrectionVectors::operator()() const
{
    if (!skew())
    {
        FatalErrorIn("skewCorrectionVectors::operator()()")
            << "Cannot return correctionVectors; mesh is not skewed"
            << abort(FatalError);
    }

    return *skewCorrectionVectors_;
}


// Do what is neccessary if the mesh has moved
bool Foam::skewCorrectionVectors::movePoints()
{
    skew_ = true;
    deleteDemandDrivenData(skewCorrectionVectors_);

    return true;
}


// ************************************************************************* //
