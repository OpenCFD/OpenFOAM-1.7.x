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

#include "extendedLeastSquaresVectors.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::extendedLeastSquaresVectors::extendedLeastSquaresVectors
(
    const fvMesh& mesh,
    const scalar minDet
)
:
    MeshObject<fvMesh, extendedLeastSquaresVectors>(mesh),
    minDet_(minDet),
    pVectorsPtr_(NULL),
    nVectorsPtr_(NULL),
    additionalCellsPtr_(NULL),
    additionalVectorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::extendedLeastSquaresVectors::~extendedLeastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    deleteDemandDrivenData(additionalCellsPtr_);
    deleteDemandDrivenData(additionalVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extendedLeastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "extendedLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    const fvMesh& mesh = mesh_;

    pVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // Build the d-vectors
    surfaceVectorField d = mesh_.Sf()/(mesh_.magSf()*mesh_.deltaCoeffs());

    if (!mesh_.orthogonal())
    {
        d -= mesh_.correctionVectors()/mesh_.deltaCoeffs();
    }


    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nCells(), symmTensor::zero);

    forAll(owner, faceI)
    {
        symmTensor wdd = 1.0/magSqr(d[faceI])*sqr(d[faceI]);

        dd[owner[faceI]] += wdd;
        dd[neighbour[faceI]] += wdd;
    }

    // Visit the boundaries. Coupled boundaries are taken into account
    // in the construction of d vectors.  
    forAll(d.boundaryField(), patchI)
    {
        const fvsPatchVectorField& pd = d.boundaryField()[patchI];

        const fvPatch& p = pd.patch();
        const unallocLabelList& faceCells = p.faceCells();

        forAll(pd, patchFaceI)
        {
            dd[faceCells[patchFaceI]] +=
                (1.0/magSqr(pd[patchFaceI]))*sqr(pd[patchFaceI]);
        }
    }

    scalarField detdd = det(dd);

    Info<< "max(detdd) = " << max(detdd) << endl;
    Info<< "min(detdd) = " << min(detdd) << endl;
    Info<< "average(detdd) = " << average(detdd) << endl;

    label nAddCells = 0;
    label maxNaddCells = 4*detdd.size();
    additionalCellsPtr_ = new List<labelPair>(maxNaddCells);
    List<labelPair>& additionalCells_ = *additionalCellsPtr_;

    forAll (detdd, i)
    {
        label count = 0;

        while (++count < 100 && detdd[i] < minDet_)
        {
            if (nAddCells == maxNaddCells)
            {
                FatalErrorIn
                (
                    "extendedLeastSquaresVectors::"
                    "makeLeastSquaresVectors() const"
                )   << "nAddCells exceeds maxNaddCells"
                    << exit(FatalError);
            }

            labelList pointLabels = mesh.cells()[i].labels(mesh.faces());

            scalar maxDetddij = 0.0;

            label addCell = -1;

            forAll(pointLabels, j)
            {
                forAll(mesh.pointCells()[pointLabels[j]], k)
                {
                    label cellj = mesh.pointCells()[pointLabels[j]][k];

                    if (cellj != i)
                    {
                        if (cellj != -1)
                        {
                            vector dCij = (mesh.C()[cellj] - mesh.C()[i]);
                            
                            symmTensor ddij = 
                                dd[i] + (1.0/magSqr(dCij))*sqr(dCij);

                            scalar detddij = det(ddij);

                            if (detddij > maxDetddij)
                            {
                                addCell = cellj;
                                maxDetddij = detddij;
                            }
                        }
                    }
                }
            }

            if (addCell != -1)
            {
                additionalCells_[nAddCells][0] = i;
                additionalCells_[nAddCells++][1] = addCell;
                vector dCij = mesh.C()[addCell] - mesh.C()[i];
                dd[i] += (1.0/magSqr(dCij))*sqr(dCij);
                detdd[i] = det(dd[i]);
            }
        }
    }

    additionalCells_.setSize(nAddCells);

    Info<< "max(detdd) = " << max(detdd) << endl;
    Info<< "min(detdd) = " << min(detdd) << endl;
    Info<< "average(detdd) = " << average(detdd) << endl;
    Info<< "nAddCells/nCells = " << scalar(nAddCells)/mesh.nCells() << endl;

    // Invert the dd tensor
    symmTensorField invDd = inv(dd);


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, faceI)
    {
        lsP[faceI] =
            (1.0/magSqr(d[faceI]))*(invDd[owner[faceI]] & d[faceI]);

        lsN[faceI] =
            ((-1.0)/magSqr(d[faceI]))*(invDd[neighbour[faceI]] & d[faceI]);
    }

    forAll(lsP.boundaryField(), patchI)
    {
        const fvsPatchVectorField& pd = d.boundaryField()[patchI];

        fvsPatchVectorField& patchLsP = lsP.boundaryField()[patchI];

        const fvPatch& p = patchLsP.patch();
        const unallocLabelList& faceCells = p.faceCells();

        forAll(p, patchFaceI)
        {
            patchLsP[patchFaceI] =
                (1.0/magSqr(pd[patchFaceI]))
               *(invDd[faceCells[patchFaceI]] & pd[patchFaceI]);
        }
    }


    additionalVectorsPtr_ = new vectorField(additionalCells_.size());
    vectorField& additionalVectors_ = *additionalVectorsPtr_;

    forAll(additionalCells_, i)
    {
        vector dCij = 
            mesh.C()[additionalCells_[i][1]] - mesh.C()[additionalCells_[i][0]];

        additionalVectors_[i] =
            (1.0/magSqr(dCij))*(invDd[additionalCells_[i][0]] & dCij);
    }

    if (debug)
    {
        Info<< "extendedLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::surfaceVectorField&
Foam::extendedLeastSquaresVectors::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField&
Foam::extendedLeastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


const Foam::List<Foam::labelPair>&
Foam::extendedLeastSquaresVectors::additionalCells() const
{
    if (!additionalCellsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *additionalCellsPtr_;
}


const Foam::vectorField&
Foam::extendedLeastSquaresVectors::additionalVectors() const
{
    if (!additionalVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *additionalVectorsPtr_;
}


bool Foam::extendedLeastSquaresVectors::movePoints()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    deleteDemandDrivenData(additionalCellsPtr_);
    deleteDemandDrivenData(additionalVectorsPtr_);

    return true;
}


// ************************************************************************* //
