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

#include "quadraticFitSnGradData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quadraticFitSnGradData, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::quadraticFitSnGradData::quadraticFitSnGradData
(
    const fvMesh& mesh,
    const scalar cWeight
)
:
    MeshObject<fvMesh, quadraticFitSnGradData>(mesh),
    centralWeight_(cWeight),
    #ifdef SPHERICAL_GEOMETRY
        dim_(2),
    #else
        dim_(mesh.nGeometricD()),
    #endif
    minSize_
    (
        dim_ == 1 ? 3 :
        dim_ == 2 ? 6 :
        dim_ == 3 ? 9 : 0
    ),
    stencil_(mesh),
    fit_(mesh.nInternalFaces())
{
    if (debug)
    {
        Info << "Contructing quadraticFitSnGradData" << endl;
    }

    // check input
    if (centralWeight_ < 1 - SMALL)
    {
        FatalErrorIn("quadraticFitSnGradData::quadraticFitSnGradData")
            << "centralWeight requested = " << centralWeight_
            << " should not be less than one"
            << exit(FatalError);
    }

    if (minSize_ == 0)
    {
        FatalErrorIn("quadraticFitSnGradData")
            << " dimension must be 1,2 or 3, not" << dim_ << exit(FatalError);
    }

    // store the polynomial size for each face to write out
    surfaceScalarField snGradPolySize
    (
        IOobject
        (
            "quadraticFitSnGradPolySize",
            "constant",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("quadraticFitSnGradPolySize", dimless, scalar(0))
    );

    // Get the cell/face centres in stencil order.
    // Centred face stencils no good for triangles of tets. Need bigger stencils
    List<List<point> > stencilPoints(stencil_.stencil().size());
    stencil_.collectData
    (
        mesh.C(),
        stencilPoints
    );

    // find the fit coefficients for every face in the mesh

    for(label faci = 0; faci < mesh.nInternalFaces(); faci++)
    {
        snGradPolySize[faci] = calcFit(stencilPoints[faci], faci);
    }

    if (debug)
    {
        snGradPolySize.write();
        Info<< "quadraticFitSnGradData::quadraticFitSnGradData() :"
            << "Finished constructing polynomialFit data"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::quadraticFitSnGradData::findFaceDirs
(
    vector& idir,        // value changed in return
    vector& jdir,        // value changed in return
    vector& kdir,        // value changed in return
    const fvMesh& mesh,
    const label faci
)
{
    idir = mesh.Sf()[faci];
    idir /= mag(idir);

    #ifndef SPHERICAL_GEOMETRY
        if (mesh.nGeometricD() <= 2) // find the normal direcion
        {
            if (mesh.geometricD()[0] == -1)
            {
                kdir = vector(1, 0, 0);
            }
            else if (mesh.geometricD()[1] == -1)
            {
                kdir = vector(0, 1, 0);
            }
            else
            {
                kdir = vector(0, 0, 1);
            }
        }
        else // 3D so find a direction in the plane of the face
        {
            const face& f = mesh.faces()[faci];
            kdir = mesh.points()[f[0]] - mesh.points()[f[1]];
        }
    #else
        // Spherical geometry so kdir is the radial direction
        kdir = mesh.Cf()[faci];
    #endif

    if (mesh.nGeometricD() == 3)
    {
        // Remove the idir component from kdir and normalise
        kdir -= (idir & kdir)*idir;

        scalar magk = mag(kdir);

        if (magk < SMALL)
        {
            FatalErrorIn("findFaceDirs") << " calculated kdir = zero"
                << exit(FatalError);
        }
        else
        {
            kdir /= magk;
        }
    }

    jdir = kdir ^ idir;
}


Foam::label Foam::quadraticFitSnGradData::calcFit
(
    const List<point>& C,
    const label faci
)
{
    vector idir(1,0,0);
    vector jdir(0,1,0);
    vector kdir(0,0,1);
    findFaceDirs(idir, jdir, kdir, mesh(), faci);

    scalarList wts(C.size(), scalar(1));
    wts[0] = centralWeight_;
    wts[1] = centralWeight_;

    point p0 = mesh().faceCentres()[faci];
    scalar scale = 0;

    // calculate the matrix of the polynomial components
    scalarRectangularMatrix B(C.size(), minSize_, scalar(0));

    for(label ip = 0; ip < C.size(); ip++)
    {
        const point& p = C[ip];

        scalar px = (p - p0)&idir;
        scalar py = (p - p0)&jdir;
        #ifdef SPHERICAL_GEOMETRY
            scalar pz = mag(p) - mag(p0);
        #else
            scalar pz = (p - p0)&kdir;
        #endif

        if (ip == 0) scale = max(max(mag(px), mag(py)), mag(pz));

        px /= scale;
        py /= scale;
        pz /= scale;

        label is = 0;

        B[ip][is++] = wts[0]*wts[ip];
        B[ip][is++] = wts[0]*wts[ip]*px;
        B[ip][is++] = wts[ip]*sqr(px);

        if (dim_ >= 2)
        {
            B[ip][is++] = wts[ip]*py;
            B[ip][is++] = wts[ip]*px*py;
            B[ip][is++] = wts[ip]*sqr(py);
        }
        if (dim_ == 3)
        {
            B[ip][is++] = wts[ip]*pz;
            B[ip][is++] = wts[ip]*px*pz;
            //B[ip][is++] = wts[ip]*py*pz;
            B[ip][is++] = wts[ip]*sqr(pz);
        }
    }

    // Set the fit
    label stencilSize = C.size();
    fit_[faci].setSize(stencilSize);
    scalarList singVals(minSize_);
    label nSVDzeros = 0;

    const scalar& deltaCoeff = mesh().deltaCoeffs()[faci];

    bool goodFit = false;
    for(int iIt = 0; iIt < 10 && !goodFit; iIt++)
    {
        SVD svd(B, SMALL);

        scalar fit0 = wts[0]*wts[0]*svd.VSinvUt()[1][0]/scale;
        scalar fit1 = wts[0]*wts[1]*svd.VSinvUt()[1][1]/scale;

        goodFit =
            fit0 < 0 && fit1 > 0
         && mag(fit0 + deltaCoeff) < 0.5*deltaCoeff
         && mag(fit1 - deltaCoeff) < 0.5*deltaCoeff;

        if (goodFit)
        {
            fit_[faci][0] = fit0;
            fit_[faci][1] = fit1;
            for(label i = 2; i < stencilSize; i++)
            {
                fit_[faci][i] = wts[0]*wts[i]*svd.VSinvUt()[1][i]/scale;
            }
            singVals = svd.S();
            nSVDzeros = svd.nZeros();
        }
        else // (not good fit so increase weight in the centre and for linear)
        {
            wts[0] *= 10;
            wts[1] *= 10;

            for(label i = 0; i < B.n(); i++)
            {
                B[i][0] *= 10;
                B[i][1] *= 10;
            }

            for(label j = 0; j < B.m(); j++)
            {
                B[0][j] *= 10;
                B[1][j] *= 10;
            }
        }
    }

    if (goodFit)
    {
        // remove the uncorrected snGradScheme coefficients
        fit_[faci][0] += deltaCoeff;
        fit_[faci][1] -= deltaCoeff;
    }
    else
    {
        Pout<< "quadratifFitSnGradData could not fit face " << faci
            << " fit_[faci][0] =  " << fit_[faci][0]
            << " fit_[faci][1] =  " << fit_[faci][1]
            << " deltaCoeff =  " << deltaCoeff << endl;
        fit_[faci] = 0;
    }

    return minSize_ - nSVDzeros;
}

bool Foam::quadraticFitSnGradData::movePoints()
{
    notImplemented("quadraticFitSnGradData::movePoints()");

    return true;
}


// ************************************************************************* //
