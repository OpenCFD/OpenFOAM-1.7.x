/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

Class
    Foam::Polynomial

Description
    Polynomial templated on size (order):

        poly = logCoeff*log(x) + sum(coeff_[i]*x^i)

    where 0 <= i <= n

    - integer powers, starting at zero
    - evaluate(x) to evaluate the poly for a given value
    - integrate(x1, x2) between two scalar values
    - integrate() to return a new, intergated coeff polynomial
      - increases the size (order)
    - integrateMinus1() to return a new, integrated coeff polynomial where
      the base poly starts at order -1

SourceFiles
    Polynomial.C

\*---------------------------------------------------------------------------*/

#ifndef Polynomial_H
#define Polynomial_H

#include "word.H"
#include "scalar.H"
#include "Ostream.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<int PolySize>
class Polynomial;

// Forward declaration of friend functions
template<int PolySize>
Ostream& operator<<
(
    Ostream&,
    const Polynomial<PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class Polynomial Declaration
\*---------------------------------------------------------------------------*/

template<int PolySize>
class Polynomial
:
    public VectorSpace<Polynomial<PolySize>, scalar, PolySize>
{
    // Private data

        //- Include the log term? - only activated using integrateMinus1()
        bool logActive_;

        //- Log coefficient - only activated using integrateMinus1()
        scalar logCoeff_;


public:

    typedef Polynomial<PolySize> polyType;

    typedef Polynomial<PolySize+1> intPolyType;


    // Constructors

        //- Construct null
        Polynomial();

        //- Construct from name and Istream
        Polynomial(const word& name, Istream& is);

        //- Copy constructor
        Polynomial(const Polynomial& poly);


    // Member Functions

        // Access

            //- Return access to the log term active flag
            bool& logActive();

            //- Return access to the log coefficient
            scalar& logCoeff();


        // Evaluation

            //- Return polynomial value
            scalar evaluate(const scalar x) const;

            //- Return integrated polynomial coefficients
            //  argument becomes zeroth element (constant of integration)
            intPolyType integrate(const scalar intConstant = 0.0);

            //- Return integrated polynomial coefficients when lowest order
            //  is -1. Argument added to zeroth element
            polyType integrateMinus1(const scalar intConstant = 0.0);

            //- Integrate between two values
            scalar integrateLimits(const scalar x1, const scalar x2) const;


    //- Ostream Operator
    friend Ostream& operator<< <PolySize>
    (
        Ostream&,
        const Polynomial&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Polynomial.C"
#   include "PolynomialIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
