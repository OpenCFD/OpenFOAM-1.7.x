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
    Foam::incompressible

Description
    Incompressible gas/liquid equation of state.

SourceFiles
    incompressibleI.H
    incompressible.C

\*---------------------------------------------------------------------------*/

#ifndef incompressible_H
#define incompressible_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class incompressible Declaration
\*---------------------------------------------------------------------------*/

class incompressible
:
    public specie
{
    // Private data

        //- Density
        scalar rho_;


public:

    // Constructors

        //- Construct from components
        inline incompressible(const specie& sp, const scalar rho);

        //- Construct from Istream
        incompressible(Istream&);

        //- Construct as named copy
        inline incompressible(const word& name, const incompressible&);

        //- Construct and return a clone
        inline autoPtr<incompressible> clone() const;

        // Selector from Istream
        inline static autoPtr<incompressible> New(Istream& is);


    // Member functions

        //- Return density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //- Return compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;

        //- Return compression factor []
        inline scalar Z(scalar p, scalar T) const;


    // Member operators

        inline void operator+=(const incompressible&);
        inline void operator-=(const incompressible&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend incompressible operator+
        (
            const incompressible&,
            const incompressible&
        );

        inline friend incompressible operator-
        (
            const incompressible&,
            const incompressible&
        );

        inline friend incompressible operator*
        (
            const scalar s,
            const incompressible&
        );

        inline friend incompressible operator==
        (
            const incompressible&,
            const incompressible&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const incompressible&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "incompressibleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
