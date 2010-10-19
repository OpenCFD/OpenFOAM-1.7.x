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
    Foam::hPolynomialThermo

Description
    Thermodynamics package templated on the equation of state, using polynomial
    functions for cp, h and s

    Polynomials for h and s derived from cp

SourceFiles
    hPolynomialThermoI.H
    hPolynomialThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hPolynomialThermo_H
#define hPolynomialThermo_H

#include "scalar.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState, int PolySize>
class hPolynomialThermo;

template<class EquationOfState, int PolySize>
inline hPolynomialThermo<EquationOfState, PolySize> operator+
(
    const hPolynomialThermo<EquationOfState, PolySize>&,
    const hPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline hPolynomialThermo<EquationOfState, PolySize> operator-
(
    const hPolynomialThermo<EquationOfState, PolySize>&,
    const hPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline hPolynomialThermo<EquationOfState, PolySize> operator*
(
    const scalar,
    const hPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline hPolynomialThermo<EquationOfState, PolySize> operator==
(
    const hPolynomialThermo<EquationOfState, PolySize>&,
    const hPolynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const hPolynomialThermo<EquationOfState, PolySize>&
);


/*---------------------------------------------------------------------------*\
                      Class hPolynomialThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState, int PolySize>
class hPolynomialThermo
:
    public EquationOfState
{
    // Private data

        //- Heat of formation
        //  Note: input in [J/kg], but internally uses [J/kmol]
        scalar Hf_;

        //- Standard entropy
        //  Note: input in [J/kg/K], but internally uses [J/kmol/K]
        scalar Sf_;

        //- Specific heat at constant pressure [J/(kg.K)]
        Polynomial<PolySize> CpPolynomial_;

        //- Enthalpy - derived from cp [J/kg] - relative to Tstd
        typename Polynomial<PolySize>::intPolyType hPolynomial_;

        //- Entropy - derived from cp [J/(kg.K)] - relative to Tstd
        Polynomial<PolySize> sPolynomial_;


    // Private member functions

        //- Construct from components
        inline hPolynomialThermo
        (
            const EquationOfState& pt,
            const scalar Hf,
            const scalar Sf,
            const Polynomial<PolySize>& cpPoly,
            const typename Polynomial<PolySize>::intPolyType& hPoly,
            const Polynomial<PolySize>& sPoly
        );


public:

    // Constructors

        //- Construct from dictionary
        hPolynomialThermo(Istream& is);

        //- Construct as copy
        inline hPolynomialThermo(const hPolynomialThermo&);

        //- Construct as a named copy
        inline hPolynomialThermo(const word&, const hPolynomialThermo&);


    // Member Functions

        //- Heat capacity at constant pressure [J/(kmol K)]
        inline scalar cp(const scalar T) const;

        //- Enthalpy [J/kmol]
        inline scalar h(const scalar T) const;

        //- Sensible enthalpy [J/kmol]
        inline scalar hs(const scalar T) const;

        //- Chemical enthalpy [J/kmol]
        inline scalar hc() const;

        //- Entropy [J/(kmol K)]
        inline scalar s(const scalar T) const;


    // Member operators

        inline hPolynomialThermo& operator=(const hPolynomialThermo&);
        inline void operator+=(const hPolynomialThermo&);
        inline void operator-=(const hPolynomialThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend hPolynomialThermo operator+ <EquationOfState, PolySize>
        (
            const hPolynomialThermo&,
            const hPolynomialThermo&
        );

        friend hPolynomialThermo operator- <EquationOfState, PolySize>
        (
            const hPolynomialThermo&,
            const hPolynomialThermo&
        );

        friend hPolynomialThermo operator* <EquationOfState, PolySize>
        (
            const scalar,
            const hPolynomialThermo&
        );

        friend hPolynomialThermo operator== <EquationOfState, PolySize>
        (
            const hPolynomialThermo&,
            const hPolynomialThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState, PolySize>
        (
            Ostream&,
            const hPolynomialThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hPolynomialThermoI.H"

#ifdef NoRepository
#   include "hPolynomialThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
