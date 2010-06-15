/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::icoPolynomial

Description
    Incompressible, polynomial form of equation of state, using a polynomial
    function for density.

SourceFiles
    icoPolynomialI.H
    icoPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef icoPolynomial_H
#define icoPolynomial_H

#include "specie.H"
#include "autoPtr.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<int PolySize>
class icoPolynomial;

template<int PolySize>
icoPolynomial<PolySize> operator+
(
    const icoPolynomial<PolySize>&,
    const icoPolynomial<PolySize>&
);

template<int PolySize>
icoPolynomial<PolySize> operator-
(
    const icoPolynomial<PolySize>&,
    const icoPolynomial<PolySize>&
);

template<int PolySize>
icoPolynomial<PolySize> operator*
(
    const scalar,
    const icoPolynomial<PolySize>&
);

template<int PolySize>
icoPolynomial<PolySize> operator==
(
    const icoPolynomial<PolySize>&,
    const icoPolynomial<PolySize>&
);

template<int PolySize>
Ostream& operator<<
(
    Ostream&,
    const icoPolynomial<PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class icoPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<int PolySize>
class icoPolynomial
:
    public specie
{
    // Private data

        //- Density
        //  Note: input in [kg/m3], but internally uses [kg/m3/kmol]
        Polynomial<PolySize> rhoPolynomial_;


public:

    // Constructors

        //- Construct from components
        inline icoPolynomial
        (
            const specie& sp,
            const Polynomial<PolySize>& rhoPoly
        );

        //- Construct from Istream
        icoPolynomial(Istream&);

        //- Construct as copy
        inline icoPolynomial(const icoPolynomial&);

        //- Construct as named copy
        inline icoPolynomial(const word& name, const icoPolynomial&);

        //- Construct and return a clone
        inline autoPtr<icoPolynomial> clone() const;

        // Selector from Istream
        inline static autoPtr<icoPolynomial> New(Istream& is);


    // Member functions

        //- Return density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;

        //- Return compression factor []
        inline scalar Z(scalar p, scalar T) const;


    // Member operators

        inline icoPolynomial& operator=(const icoPolynomial&);
        inline void operator+=(const icoPolynomial&);
        inline void operator-=(const icoPolynomial&);

        inline void operator*=(const scalar);


    // Friend operators

        friend icoPolynomial operator+ <PolySize>
        (
            const icoPolynomial&,
            const icoPolynomial&
        );

        friend icoPolynomial operator- <PolySize>
        (
            const icoPolynomial&,
            const icoPolynomial&
        );

        friend icoPolynomial operator* <PolySize>
        (
            const scalar s,
            const icoPolynomial&
        );

        friend icoPolynomial operator== <PolySize>
        (
            const icoPolynomial&,
            const icoPolynomial&
        );


    // Ostream Operator

        friend Ostream& operator<< <PolySize>(Ostream&, const icoPolynomial&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIcoPolynomial(PolySize)                                          \
                                                                             \
defineTemplateTypeNameAndDebugWithName                                       \
(                                                                            \
    icoPolynomial<PolySize>,                                                 \
    "icoPolynomial<"#PolySize">",                                            \
    0                                                                        \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "icoPolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "icoPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
