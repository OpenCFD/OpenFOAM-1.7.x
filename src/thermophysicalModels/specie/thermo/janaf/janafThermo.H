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

Class
    Foam::janafThermo

Description
    JANAF tables based thermodynamics package templated
    into the equationOfState.

SourceFiles
    janafThermoI.H
    janafThermo.C

\*---------------------------------------------------------------------------*/

#ifndef janafThermo_H
#define janafThermo_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class equationOfState> class janafThermo;

template<class equationOfState>
inline janafThermo<equationOfState> operator+
(
    const janafThermo<equationOfState>&,
    const janafThermo<equationOfState>&
);

template<class equationOfState>
inline janafThermo<equationOfState> operator-
(
    const janafThermo<equationOfState>&,
    const janafThermo<equationOfState>&
);

template<class equationOfState>
inline janafThermo<equationOfState> operator*
(
    const scalar,
    const janafThermo<equationOfState>&
);

template<class equationOfState>
inline janafThermo<equationOfState> operator==
(
    const janafThermo<equationOfState>&,
    const janafThermo<equationOfState>&
);

template<class equationOfState>
Ostream& operator<<
(
    Ostream&,
    const janafThermo<equationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class janafThermo Declaration
\*---------------------------------------------------------------------------*/

template<class equationOfState>
class janafThermo
:
    public equationOfState
{

public:

        static const int nCoeffs_ = 7;
        typedef scalar coeffArray[7];

private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;


    // Private member functions

        //- Check given temperature is within the range of the fitted coeffs
        inline void checkT(const scalar T) const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline janafThermo
        (
            const equationOfState& st,
            const scalar Tlow,
            const scalar Thigh,
            const scalar Tcommon,
            const coeffArray& highCpCoeffs,
            const coeffArray& lowCpCoeffs
        );

        //- Construct from Istream
        janafThermo(Istream&);

        //- Construct as a named copy
        inline janafThermo(const word&, const janafThermo&);


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

        inline void operator+=(const janafThermo&);
        inline void operator-=(const janafThermo&);


    // Friend operators

        friend janafThermo operator+ <equationOfState>
        (
            const janafThermo&,
            const janafThermo&
        );

        friend janafThermo operator- <equationOfState>
        (
            const janafThermo&,
            const janafThermo&
        );

        friend janafThermo operator* <equationOfState>
        (
            const scalar,
            const janafThermo&
        );

        friend janafThermo operator== <equationOfState>
        (
            const janafThermo&,
            const janafThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <equationOfState>
        (
            Ostream&,
            const janafThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "janafThermoI.H"

#ifdef NoRepository
#   include "janafThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
