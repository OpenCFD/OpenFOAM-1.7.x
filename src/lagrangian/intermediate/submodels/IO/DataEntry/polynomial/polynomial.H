/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::polynomial

Description
    Polynomial container data entry for scalars. Items are stored in a list of
    Tuple2's. Data is input in the form, e.g. for an entry <entryName> that
    describes y = x^2 + 2x^3

    @verbatim
        <entryName>   polynomial
        (
            (1    2)
            (2    3)
        );
    @endverbatim

SourceFiles
    polynomial.C

\*---------------------------------------------------------------------------*/

#ifndef polynomial_H
#define polynomial_H

#include "DataEntry.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polynomial;

// Forward declaration of friend functions
Ostream& operator<<
(
    Ostream&,
    const polynomial&
);

/*---------------------------------------------------------------------------*\
                        Class polynomial Declaration
\*---------------------------------------------------------------------------*/

class polynomial
:
    public DataEntry<scalar>
{
    // Private data

        //- Polynomial coefficients - list of prefactor, exponent
        List<Tuple2<scalar, scalar> > coeffs_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const polynomial&);


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from entry name and Istream
        polynomial(const word& entryName, Istream& is);

        //- Copy constructor
        polynomial(const polynomial& poly);


    //- Destructor
    virtual ~polynomial();


    // Member Functions

        //- Return polynomial value
        scalar value(const scalar x) const;

        //- Integrate between two (scalar) values
        scalar integrate(const scalar x1, const scalar x2) const;


    //- Ostream Operator
    friend Ostream& operator<<
    (
        Ostream&,
        const polynomial&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
