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

#include "GaussSeidelSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GaussSeidelSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<GaussSeidelSmoother>
        addGaussSeidelSmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable<GaussSeidelSmoother>
        addGaussSeidelSmootherAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussSeidelSmoother::GaussSeidelSmoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GaussSeidelSmoother::smooth
(
    const word& fieldName_,
    scalarField& psi,
    const lduMatrix& matrix_,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs_,
    const lduInterfaceFieldPtrsList& interfaces_,
    const direction cmpt,
    const label nSweeps
)
{
    register scalar* __restrict__ psiPtr = psi.begin();

    register const label nCells = psi.size();

    scalarField bPrime(nCells);
    register scalar* __restrict__ bPrimePtr = bPrime.begin();

    register const scalar* const __restrict__ diagPtr = matrix_.diag().begin();
    register const scalar* const __restrict__ upperPtr =
        matrix_.upper().begin();
    register const scalar* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    register const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    register const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();


    // Parallel boundary initialisation.  The parallel boundary is treated
    // as an effective jacobi interface in the boundary.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    FieldField<Field, scalar> mBouCoeffs(interfaceBouCoeffs_.size());

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs.set(patchi, -interfaceBouCoeffs_[patchi]);
        }
    }

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        bPrime = source;

        matrix_.initMatrixInterfaces
        (
            mBouCoeffs,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );

        matrix_.updateMatrixInterfaces
        (
            mBouCoeffs,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );

        register scalar curPsi;
        register label fStart;
        register label fEnd = ownStartPtr[0];

        for (register label cellI=0; cellI<nCells; cellI++)
        {
            // Start and end of this row
            fStart = fEnd;
            fEnd = ownStartPtr[cellI + 1];

            // Get the accumulated neighbour side
            curPsi = bPrimePtr[cellI];

            // Accumulate the owner product side
            for (register label curFace=fStart; curFace<fEnd; curFace++)
            {
                curPsi -= upperPtr[curFace]*psiPtr[uPtr[curFace]];
            }

            // Finish current psi
            curPsi /= diagPtr[cellI];

            // Distribute the neighbour side using current psi
            for (register label curFace=fStart; curFace<fEnd; curFace++)
            {
                bPrimePtr[uPtr[curFace]] -= lowerPtr[curFace]*curPsi;
            }

            psiPtr[cellI] = curPsi;
        }
    }
}


void Foam::GaussSeidelSmoother::smooth
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        fieldName_,
        psi,
        matrix_,
        source,
        interfaceBouCoeffs_,
        interfaces_,
        cmpt,
        nSweeps
    );
}


// ************************************************************************* //
