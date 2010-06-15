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
    Foam::dimensionSet

Description
    Dimension set for the base types.
    This type may be used to implement rigorous dimension checking
    for algebraic manipulation.

SourceFiles
    dimensionSet.C
    dimensionSetIO.C
    dimensionSets.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionSet_H
#define dimensionSet_H

#include "scalar.H"
#include "bool.H"
#include "dimensionedScalarFwd.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class dimensionSet;

// Friend functions

dimensionSet max(const dimensionSet&, const dimensionSet&);
dimensionSet min(const dimensionSet&, const dimensionSet&);
dimensionSet cmptMultiply(const dimensionSet&, const dimensionSet&);
dimensionSet cmptDivide(const dimensionSet&, const dimensionSet&);

dimensionSet pow(const dimensionSet&, const scalar);
dimensionSet pow(const dimensionSet&, const dimensionedScalar&);
dimensionSet pow(const dimensionedScalar&, const dimensionSet&);

dimensionSet sqr(const dimensionSet&);
dimensionSet pow3(const dimensionSet&);
dimensionSet pow4(const dimensionSet&);
dimensionSet pow5(const dimensionSet&);
dimensionSet pow6(const dimensionSet&);

dimensionSet sqrt(const dimensionSet&);
dimensionSet magSqr(const dimensionSet&);
dimensionSet mag(const dimensionSet&);
dimensionSet sign(const dimensionSet&);
dimensionSet pos(const dimensionSet&);
dimensionSet neg(const dimensionSet&);
dimensionSet inv(const dimensionSet&);

// Function to check the argument is dimensionless
//  for transcendental functions
dimensionSet trans(const dimensionSet&);

// Return the argument; transformations do not change the dimensions
dimensionSet transform(const dimensionSet&);

// Friend operators

dimensionSet operator-(const dimensionSet&);
dimensionSet operator+(const dimensionSet&, const dimensionSet&);
dimensionSet operator-(const dimensionSet&, const dimensionSet&);
dimensionSet operator*(const dimensionSet&, const dimensionSet&);
dimensionSet operator/(const dimensionSet&, const dimensionSet&);
dimensionSet operator&(const dimensionSet&, const dimensionSet&);
dimensionSet operator^(const dimensionSet&, const dimensionSet&);
dimensionSet operator&&(const dimensionSet&, const dimensionSet&);

// IOstream operators

Istream& operator>>(Istream&, dimensionSet&);
Ostream& operator<<(Ostream&, const dimensionSet&);


/*---------------------------------------------------------------------------*\
                           Class dimensionSet Declaration
\*---------------------------------------------------------------------------*/

class dimensionSet
{

public:

    // Member constants

        enum
        {
            nDimensions = 7    // Number of dimensions in SI is 7
        };

        //- Define an enumeration for the names of the dimension exponents
        enum dimensionType
        {
            MASS,               // kilogram   kg
            LENGTH,             // metre      m
            TIME,               // second     s
            TEMPERATURE,        // Kelvin     K
            MOLES,              // mole       mol
            CURRENT,            // Ampere     A
            LUMINOUS_INTENSITY  // Candela    Cd
        };


    // Static data members

        static const scalar smallExponent;


private:

    // private data

        // dimensionSet stored as an array of dimension exponents
        scalar exponents_[nDimensions];


public:

    // Declare name of the class and its debug switch
    ClassName("dimensionSet");


    // Constructors

        //- Construct given individual dimension exponents for all
        //  seven dimensions
        dimensionSet
        (
            const scalar mass,
            const scalar length,
            const scalar time,
            const scalar temperature,
            const scalar moles,
            const scalar current,
            const scalar luminousIntensity
        );

        //- Construct given individual dimension exponents for first
        //  five dimensions
        dimensionSet
        (
            const scalar mass,
            const scalar length,
            const scalar time,
            const scalar temperature,
            const scalar moles
        );

        //- Construct from Istream
        dimensionSet(Istream&);


    // Member functions

        bool dimensionless() const;
        void reset(const dimensionSet&);


    // Member operators

        scalar operator[](const dimensionType) const;
        scalar& operator[](const dimensionType);
        bool operator==(const dimensionSet&) const;
        bool operator!=(const dimensionSet&) const;

        bool operator=(const dimensionSet&) const;

        bool operator+=(const dimensionSet&) const;
        bool operator-=(const dimensionSet&) const;
        bool operator*=(const dimensionSet&);
        bool operator/=(const dimensionSet&);


    // Friend functions

        friend dimensionSet max(const dimensionSet&, const dimensionSet&);
        friend dimensionSet min(const dimensionSet&, const dimensionSet&);
        friend dimensionSet cmptMultiply
        (
            const dimensionSet&,
            const dimensionSet&
        );
        friend dimensionSet cmptDivide
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet pow(const dimensionSet&, const scalar);
        friend dimensionSet pow(const dimensionSet&, const dimensionedScalar&);
        friend dimensionSet pow(const dimensionedScalar&, const dimensionSet&);

        friend dimensionSet sqr(const dimensionSet&);
        friend dimensionSet pow3(const dimensionSet&);
        friend dimensionSet pow4(const dimensionSet&);
        friend dimensionSet pow5(const dimensionSet&);
        friend dimensionSet pow6(const dimensionSet&);

        friend dimensionSet sqrt(const dimensionSet&);
        friend dimensionSet magSqr(const dimensionSet&);
        friend dimensionSet mag(const dimensionSet&);
        friend dimensionSet sign(const dimensionSet&);
        friend dimensionSet pos(const dimensionSet&);
        friend dimensionSet neg(const dimensionSet&);
        friend dimensionSet inv(const dimensionSet&);

        //- Function to check the argument is dimensionless
        //  for transcendental functions
        friend dimensionSet trans(const dimensionSet&);

        //- Return the argument; transformations do not change the dimensions
        friend dimensionSet transform(const dimensionSet&);


    // Friend operators

        friend dimensionSet operator-(const dimensionSet&);

        friend dimensionSet operator+
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator-
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator*
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator/
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator&
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator^
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator&&
        (
            const dimensionSet&,
            const dimensionSet&
        );


    // IOstream operators

        friend Istream& operator>>(Istream&, dimensionSet&);
        friend Ostream& operator<<(Ostream&, const dimensionSet&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dimensionSets.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
