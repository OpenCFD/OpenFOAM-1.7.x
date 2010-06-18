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

Class
    Foam::fvMatrix

Description
    Finite-Volume matrix.

SourceFiles
    fvMatrix.C
    fvMatrixSolve.C
    fvScalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef fvMatrix_H
#define fvMatrix_H

#include "volFields.H"
#include "surfaceFields.H"
#include "lduMatrix.H"
#include "tmp.H"
#include "autoPtr.H"
#include "dimensionedTypes.H"
#include "zeroField.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class fvMatrix;

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const fvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const fvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<fvMatrix<Type> >&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<fvMatrix<Type> >&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<fvMatrix<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
Ostream& operator<<(Ostream&, const fvMatrix<Type>&);


/*---------------------------------------------------------------------------*\
                           Class fvMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvMatrix
:
    public refCount,
    public lduMatrix
{
    // Private data

        // Reference to GeometricField<Type, fvPatchField, volMesh>
        GeometricField<Type, fvPatchField, volMesh>& psi_;

        //- Dimension set
        dimensionSet dimensions_;

        //- Source term
        Field<Type> source_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //  for internal cells
        FieldField<Field, Type> internalCoeffs_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //  for boundary cells
        FieldField<Field, Type> boundaryCoeffs_;


        //- Face flux field for non-orthogonal correction
        mutable GeometricField<Type, fvsPatchField, surfaceMesh>
            *faceFluxCorrectionPtr_;


    // Private member functions

        //- Add patch contribution to internal field
        template<class Type2>
        void addToInternalField
        (
            const unallocLabelList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void addToInternalField
        (
            const unallocLabelList& addr,
            const tmp<Field<Type2> >& tpf,
            Field<Type2>& intf
        ) const;

        //- Subtract patch contribution from internal field
        template<class Type2>
        void subtractFromInternalField
        (
            const unallocLabelList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void subtractFromInternalField
        (
            const unallocLabelList& addr,
            const tmp<Field<Type2> >& tpf,
            Field<Type2>& intf
        ) const;


        // Matrix completion functionality

            void addBoundaryDiag
            (
                scalarField& diag,
                const direction cmpt
            ) const;

            void addCmptAvBoundaryDiag(scalarField& diag) const;

            void addBoundarySource
            (
                Field<Type>& source,
                const bool couples=true
            ) const;


public:

    //- Solver class returned by the solver function
    //  used for systems in which it is useful to cache the solver for reuse
    //  e.g. if the solver is potentialy expensive to construct (AMG) and can
    //  be used several times (PISO)
    class fvSolver
    {
        fvMatrix<Type>& fvMat_;

        autoPtr<lduMatrix::solver> solver_;

    public:

        // Constructors

            fvSolver(fvMatrix<Type>& fvMat, autoPtr<lduMatrix::solver> sol)
            :
                fvMat_(fvMat),
                solver_(sol)
            {}


        // Member functions

            //- Solve returning the solution statistics.
            //  Use the given solver controls
            lduMatrix::solverPerformance solve(const dictionary&);

            //- Solve returning the solution statistics.
            //  Solver controls read from fvSolution
            lduMatrix::solverPerformance solve();
    };


    ClassName("fvMatrix");


    // Constructors

        //- Construct given a field to solve for
        fvMatrix
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const dimensionSet&
        );

        //- Construct as copy
        fvMatrix(const fvMatrix<Type>&);

        //- Construct as copy of tmp<fvMatrix<Type> > deleting argument
#       ifdef ConstructFromTmp
        fvMatrix(const tmp<fvMatrix<Type> >&);
#       endif

        //- Construct from Istream given field to solve for
        fvMatrix(GeometricField<Type, fvPatchField, volMesh>&, Istream&);


    // Destructor

        virtual ~fvMatrix();


    // Member functions

        // Access

            const GeometricField<Type, fvPatchField, volMesh>& psi() const
            {
                return psi_;
            }

            GeometricField<Type, fvPatchField, volMesh>& psi()
            {
                return psi_;
            }

            const dimensionSet& dimensions() const
            {
                return dimensions_;
            }

            Field<Type>& source()
            {
                return source_;
            }

            const Field<Type>& source() const
            {
                return source_;
            }

            //- fvBoundary scalar field containing pseudo-matrix coeffs
            //  for internal cells
            FieldField<Field, Type>& internalCoeffs()
            {
                return internalCoeffs_;
            }

            //- fvBoundary scalar field containing pseudo-matrix coeffs
            //  for boundary cells
            FieldField<Field, Type>& boundaryCoeffs()
            {
                return boundaryCoeffs_;
            }


            //- Declare return type of the faceFluxCorrectionPtr() function
            typedef GeometricField<Type, fvsPatchField, surfaceMesh>
                *surfaceTypeFieldPtr;

            //- Return pointer to face-flux non-orthogonal correction field
            surfaceTypeFieldPtr& faceFluxCorrectionPtr()
            {
                return faceFluxCorrectionPtr_;
            }


        // Operations

            //- Set solution in given cells and eliminate corresponding
            //  equations from the matrix
            void setValues
            (
                const labelList& cells,
                const Field<Type>& values
            );

            //- Set reference level for solution
            void setReference
            (
                const label celli,
                const Type& value,
                const bool forceReference = false
            );

            //- Set reference level for a component of the solution
            //  on a given patch face
            void setComponentReference
            (
                const label patchi,
                const label facei,
                const direction cmpt,
                const scalar value
            );

            //- Relax matrix (for steady-state solution).
            //  alpha = 1 : diagonally equal
            //  alpha < 1 : diagonally dominant
            //  alpha = 0 : do nothing
            void relax(const scalar alpha);

            //- Relax matrix (for steady-state solution).
            //  alpha is read from controlDict
            void relax();

            //- Manipulate based on a boundary field
            void boundaryManipulate
            (
                typename GeometricField<Type, fvPatchField, volMesh>::
                    GeometricBoundaryField& values
            );

            //- Construct and return the solver
            //  Use the given solver controls
            autoPtr<fvSolver> solver(const dictionary&);

            //- Construct and return the solver
            //  Solver controls read from fvSolution
            autoPtr<fvSolver> solver();

            //- Solve returning the solution statistics.
            //  Use the given solver controls
            lduMatrix::solverPerformance solve(const dictionary&);

            //- Solve returning the solution statistics.
            //  Solver controls read from fvSolution
            lduMatrix::solverPerformance solve();

            //- Return the matrix residual
            tmp<Field<Type> > residual() const;

            //- Return the matrix scalar diagonal
            tmp<scalarField> D() const;

            //- Return the matrix Type diagonal
            tmp<Field<Type> > DD() const;

            //- Return the central coefficient
            tmp<volScalarField> A() const;

            //- Return the H operation source
            tmp<GeometricField<Type, fvPatchField, volMesh> > H() const;

            //- Return H(1)
            tmp<volScalarField> H1() const;

            //- Return the face-flux field from the matrix
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
                flux() const;


    // Member operators

        void operator=(const fvMatrix<Type>&);
        void operator=(const tmp<fvMatrix<Type> >&);

        void negate();

        void operator+=(const fvMatrix<Type>&);
        void operator+=(const tmp<fvMatrix<Type> >&);

        void operator-=(const fvMatrix<Type>&);
        void operator-=(const tmp<fvMatrix<Type> >&);

        void operator+=
        (
            const DimensionedField<Type, volMesh>&
        );
        void operator+=
        (
            const tmp<DimensionedField<Type, volMesh> >&
        );
        void operator+=
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        );

        void operator-=
        (
            const DimensionedField<Type, volMesh>&
        );
        void operator-=
        (
            const tmp<DimensionedField<Type, volMesh> >&
        );
        void operator-=
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        );

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator+=(const zeroField&);
        void operator-=(const zeroField&);

        void operator*=(const DimensionedField<scalar, volMesh>&);
        void operator*=(const tmp<DimensionedField<scalar, volMesh> >&);
        void operator*=(const tmp<volScalarField>&);

        void operator*=(const dimensioned<scalar>&);


    // Friend operators

        friend tmp<GeometricField<Type, fvPatchField, volMesh> >
        operator& <Type>
        (
            const fvMatrix<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        friend tmp<GeometricField<Type, fvPatchField, volMesh> >
        operator& <Type>
        (
            const fvMatrix<Type>&,
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        );

        friend tmp<GeometricField<Type, fvPatchField, volMesh> >
        operator& <Type>
        (
            const tmp<fvMatrix<Type> >&,
            const DimensionedField<Type, volMesh>&
        );

        friend tmp<GeometricField<Type, fvPatchField, volMesh> >
        operator& <Type>
        (
            const tmp<fvMatrix<Type> >&,
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        );


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const fvMatrix<Type>&
        );
};


// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const fvMatrix<Type>&,
    const fvMatrix<Type>&,
    const char*
);

template<class Type>
void checkMethod
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&,
    const char*
);

template<class Type>
void checkMethod
(
    const fvMatrix<Type>&,
    const dimensioned<Type>&,
    const char*
);


//- Solve returning the solution statistics given convergence tolerance
//  Use the given solver controls
template<class Type>
lduMatrix::solverPerformance solve(fvMatrix<Type>&, const dictionary&);


//- Solve returning the solution statistics given convergence tolerance,
//  deleting temporary matrix after solution.
//  Use the given solver controls
template<class Type>
lduMatrix::solverPerformance solve
(
    const tmp<fvMatrix<Type> >&,
    const dictionary&
);


//- Solve returning the solution statistics given convergence tolerance
//  Solver controls read fvSolution
template<class Type>
lduMatrix::solverPerformance solve(fvMatrix<Type>&);


//- Solve returning the solution statistics given convergence tolerance,
//  deleting temporary matrix after solution.
//  Solver controls read fvSolution
template<class Type>
lduMatrix::solverPerformance solve(const tmp<fvMatrix<Type> >&);


//- Return the correction form of the given matrix
//  by subtracting the matrix multiplied by the current field
template<class Type>
tmp<fvMatrix<Type> > correction(const fvMatrix<Type>&);


//- Return the correction form of the given temporary matrix
//  by subtracting the matrix multiplied by the current field
template<class Type>
tmp<fvMatrix<Type> > correction(const tmp<fvMatrix<Type> >&);


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >&,
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >&,
    const dimensioned<Type>&
);


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>&,
    const zeroField&
);

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >&,
    const zeroField&
);


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >&,
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const DimensionedField<Type, volMesh>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<DimensionedField<Type, volMesh> >&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const DimensionedField<Type, volMesh>&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<DimensionedField<Type, volMesh> >&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const dimensioned<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const dimensioned<Type>&,
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >&,
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const DimensionedField<Type, volMesh>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<DimensionedField<Type, volMesh> >&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const DimensionedField<Type, volMesh>&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<DimensionedField<Type, volMesh> >&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const dimensioned<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const dimensioned<Type>&,
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const DimensionedField<scalar, volMesh>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const tmp<DimensionedField<scalar, volMesh> >&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const tmp<volScalarField>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const DimensionedField<scalar, volMesh>&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const tmp<DimensionedField<scalar, volMesh> >&,
    const tmp<fvMatrix<Type> >&
);

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const tmp<volScalarField>&,
    const tmp<fvMatrix<Type> >&
);


template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const dimensioned<scalar>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const dimensioned<scalar>&,
    const tmp<fvMatrix<Type> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvMatrix.C"
#endif

// Specialisation for scalars
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
