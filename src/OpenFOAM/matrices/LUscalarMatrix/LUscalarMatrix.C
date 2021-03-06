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

#include "LUscalarMatrix.H"
#include "lduMatrix.H"
#include "procLduMatrix.H"
#include "procLduInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LUscalarMatrix::LUscalarMatrix(const scalarSquareMatrix& matrix)
:
    scalarSquareMatrix(matrix),
    pivotIndices_(n())
{
    LUDecompose(*this, pivotIndices_);
}


Foam::LUscalarMatrix::LUscalarMatrix
(
    const lduMatrix& ldum,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
{
    if (Pstream::parRun())
    {
        PtrList<procLduMatrix> lduMatrices(Pstream::nProcs());

        label lduMatrixi = 0;

        lduMatrices.set
        (
            lduMatrixi++,
            new procLduMatrix
            (
                ldum,
                interfaceCoeffs,
                interfaces
            )
        );

        if (Pstream::master())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                lduMatrices.set
                (
                    lduMatrixi++,
                    new procLduMatrix(IPstream(Pstream::scheduled, slave)())
                );
            }
        }
        else
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            procLduMatrix cldum
            (
                ldum,
                interfaceCoeffs,
                interfaces
            );
            toMaster<< cldum;

        }

        if (Pstream::master())
        {
            label nCells = 0;
            forAll(lduMatrices, i)
            {
                nCells += lduMatrices[i].size();
            }

            scalarSquareMatrix m(nCells, nCells, 0.0);
            transfer(m);
            convert(lduMatrices);
        }
    }
    else
    {
        label nCells = ldum.lduAddr().size();
        scalarSquareMatrix m(nCells, nCells, 0.0);
        transfer(m);
        convert(ldum, interfaceCoeffs, interfaces);
    }

    if (Pstream::master())
    {
        pivotIndices_.setSize(n());
        LUDecompose(*this, pivotIndices_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LUscalarMatrix::convert
(
    const lduMatrix& ldum,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
{
    const label* __restrict__ uPtr = ldum.lduAddr().upperAddr().begin();
    const label* __restrict__ lPtr = ldum.lduAddr().lowerAddr().begin();

    const scalar* __restrict__ diagPtr = ldum.diag().begin();
    const scalar* __restrict__ upperPtr = ldum.upper().begin();
    const scalar* __restrict__ lowerPtr = ldum.lower().begin();

    register const label nCells = ldum.diag().size();
    register const label nFaces = ldum.upper().size();

    for (register label cell=0; cell<nCells; cell++)
    {
        operator[](cell)[cell] = diagPtr[cell];
    }

    for (register label face=0; face<nFaces; face++)
    {
        label uCell = uPtr[face];
        label lCell = lPtr[face];

        operator[](uCell)[lCell] = lowerPtr[face];
        operator[](lCell)[uCell] = upperPtr[face];
    }

    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            const lduInterface& interface = interfaces[inti].interface();

            const label* __restrict__ ulPtr = interface.faceCells().begin();
            const scalar* __restrict__ upperLowerPtr =
                interfaceCoeffs[inti].begin();

            register label inFaces = interface.faceCells().size()/2;

            for (register label face=0; face<inFaces; face++)
            {
                label uCell = ulPtr[face];
                label lCell = ulPtr[face + inFaces];

                operator[](uCell)[lCell] -= upperLowerPtr[face + inFaces];
                operator[](lCell)[uCell] -= upperLowerPtr[face];
            }
        }
    }

    //printDiagonalDominance();
}


void Foam::LUscalarMatrix::convert
(
    const PtrList<procLduMatrix>& lduMatrices
)
{
    procOffsets_.setSize(lduMatrices.size() + 1);
    procOffsets_[0] = 0;

    forAll(lduMatrices, ldumi)
    {
        procOffsets_[ldumi+1] = procOffsets_[ldumi] + lduMatrices[ldumi].size();
    }

    forAll(lduMatrices, ldumi)
    {
        const procLduMatrix& lduMatrixi = lduMatrices[ldumi];
        label offset = procOffsets_[ldumi];

        const label* __restrict__ uPtr = lduMatrixi.upperAddr_.begin();
        const label* __restrict__ lPtr = lduMatrixi.lowerAddr_.begin();

        const scalar* __restrict__ diagPtr = lduMatrixi.diag_.begin();
        const scalar* __restrict__ upperPtr = lduMatrixi.upper_.begin();
        const scalar* __restrict__ lowerPtr = lduMatrixi.lower_.begin();

        register const label nCells = lduMatrixi.size();
        register const label nFaces = lduMatrixi.upper_.size();

        for (register label cell=0; cell<nCells; cell++)
        {
            label globalCell = cell + offset;
            operator[](globalCell)[globalCell] = diagPtr[cell];
        }

        for (register label face=0; face<nFaces; face++)
        {
            label uCell = uPtr[face] + offset;
            label lCell = lPtr[face] + offset;

            operator[](uCell)[lCell] = lowerPtr[face];
            operator[](lCell)[uCell] = upperPtr[face];
        }

        const PtrList<procLduInterface>& interfaces =
            lduMatrixi.interfaces_;

        forAll(interfaces, inti)
        {
            const procLduInterface& interface = interfaces[inti];

            if (interface.myProcNo_ == interface.neighbProcNo_)
            {
                const label* __restrict__ ulPtr = interface.faceCells_.begin();

                const scalar* __restrict__ upperLowerPtr =
                    interface.coeffs_.begin();

                register label inFaces = interface.faceCells_.size()/2;

                for (register label face=0; face<inFaces; face++)
                {
                    label uCell = ulPtr[face] + offset;
                    label lCell = ulPtr[face + inFaces] + offset;

                    operator[](uCell)[lCell] -= upperLowerPtr[face + inFaces];
                    operator[](lCell)[uCell] -= upperLowerPtr[face];
                }
            }
            else if (interface.myProcNo_ < interface.neighbProcNo_)
            {
                const PtrList<procLduInterface>& neiInterfaces =
                    lduMatrices[interface.neighbProcNo_].interfaces_;

                label neiInterfacei = -1;

                forAll(neiInterfaces, ninti)
                {
                    if
                    (
                        neiInterfaces[ninti].neighbProcNo_
                     == interface.myProcNo_
                    )
                    {
                        neiInterfacei = ninti;
                        break;
                    }
                }

                if (neiInterfacei == -1)
                {
                    FatalErrorIn("LUscalarMatrix::convert") << exit(FatalError);
                }

                const procLduInterface& neiInterface =
                    neiInterfaces[neiInterfacei];

                const label* __restrict__ uPtr = interface.faceCells_.begin();
                const label* __restrict__ lPtr = neiInterface.faceCells_.begin();

                const scalar* __restrict__ upperPtr = interface.coeffs_.begin();
                const scalar* __restrict__ lowerPtr = neiInterface.coeffs_.begin();

                register label inFaces = interface.faceCells_.size();
                label neiOffset = procOffsets_[interface.neighbProcNo_];

                for (register label face=0; face<inFaces; face++)
                {
                    label uCell = uPtr[face] + offset;
                    label lCell = lPtr[face] + neiOffset;

                    operator[](uCell)[lCell] -= lowerPtr[face];
                    operator[](lCell)[uCell] -= upperPtr[face];
                }
            }
        }
    }

    //printDiagonalDominance();
}


void Foam::LUscalarMatrix::printDiagonalDominance() const
{
    for (label i=0; i<n(); i++)
    {
        scalar sum = 0.0;
        for (label j=0; j<n(); j++)
        {
            if (i != j)
            {
                sum += operator[](i)[j];
            }
        }
        Info<< mag(sum)/mag(operator[](i)[i]) << endl;
    }
}


// ************************************************************************* //
