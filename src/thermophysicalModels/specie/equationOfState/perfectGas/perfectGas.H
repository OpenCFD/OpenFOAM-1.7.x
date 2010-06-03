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
    Foam::perfectGas

Description
    Perfect gas equation of state.

SourceFiles
    perfectGasI.H
    perfectGas.C

\*---------------------------------------------------------------------------*/

#ifndef perfectGas_H
#define perfectGas_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class perfectGas Declaration
\*---------------------------------------------------------------------------*/

class perfectGas
:
    public specie
{

public:

    // Constructors

        //- Construct from components
        inline perfectGas(const specie& sp);

        //- Construct from Istream
        perfectGas(Istream&);

        //- Construct as named copy
        inline perfectGas(const word& name, const perfectGas&);

        //- Construct and return a clone
        inline autoPtr<perfectGas> clone() const;

        // Selector from Istream
        inline static autoPtr<perfectGas> New(Istream& is);


    // Member functions

        //- Return density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;

        //- Return compression factor []
        inline scalar Z(scalar p, scalar T) const;


    // Member operators

        inline void operator+=(const perfectGas&);
        inline void operator-=(const perfectGas&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend perfectGas operator+
        (
            const perfectGas&,
            const perfectGas&
        );

        inline friend perfectGas operator-
        (
            const perfectGas&,
            const perfectGas&
        );

        inline friend perfectGas operator*
        (
            const scalar s,
            const perfectGas&
        );

        inline friend perfectGas operator==
        (
            const perfectGas&,
            const perfectGas&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const perfectGas&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "perfectGasI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
