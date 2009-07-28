/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::Tuple

Description
    A 2 Tuple. Differs from Tuple in that the two elements can be different
    type.


\*---------------------------------------------------------------------------*/

#ifndef Tuple_H
#define Tuple_H

#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type1, class Type2>
class Tuple;

template<class Type1, class Type2>
Istream& operator>>(Istream&, Tuple<Type1, Type2>&);

template<class Type1, class Type2>
Ostream& operator<<(Ostream&, const Tuple<Type1, Type2>&);

/*---------------------------------------------------------------------------*\
                            Class Tuple Declaration
\*---------------------------------------------------------------------------*/

template<class Type1, class Type2>
class Tuple
{
    // Private data

       Type1 first_;
       Type2 second_;


public:

    // Constructors

        //- Null constructor for lists
        inline Tuple()
        {}

        //- Construct from components
        inline Tuple(const Type1& first, const Type2& second)
        :
            first_(first),
            second_(second)
        {}

        //- Construct from Istream
        inline Tuple(Istream& is)
        {
            // Read beginning of pair
            is.readBegin("pair");

            is >> first_ >> second_;

            // Read end of pair
            is.readEnd("pair");

            // Check state of Istream
            is.check("Tuple::Tuple(Istream&)");
        }


    // Member Functions

        //- Return first
        inline Type1 first() const
        {
            return first_;
        }

        //- Return first
        inline Type1& first()
        {
            return first_;
        }

        //- Return second
        inline Type2 second() const
        {
            return second_;
        }

        //- Return second
        inline Type2& second()
        {
            return second_;
        }

        //- Return reverse pair
        inline Tuple<Type1, Type2> reverseTuple() const
        {
            return Tuple<Type1, Type2>(second_, first_);
        }


    // Friend Operators

        inline friend bool operator==
        (
            const Tuple<Type1, Type2>& a,
            const Tuple<Type1, Type2>& b
        )
        {
            return
            (
                (a.first_ == b.first_) && (a.second_ == b.second_)
            );
        }

        inline friend bool operator!=
        (
            const Tuple<Type1, Type2>& a,
            const Tuple<Type1, Type2>& b
        )
        {
            return (!(a == b));
        }


    // IOstream Operators

        friend Istream& operator>> <Type1, Type2>
        (
            Istream& is,
            Tuple<Type1, Type2>& p
        );
        friend Ostream& operator<< <Type1, Type2>
        (
            Ostream& os,
            const Tuple<Type1, Type2>& p
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type1, class Type2>
Istream& operator>>(Istream& is, Tuple<Type1, Type2>& p)
{
    // Read beginning of Tuple<Type, Type>
    is.readBegin("Tuple<Type, Type>");

    is >> p.first_ >> p.second_;

    // Read end of Tuple<Type, Type>
    is.readEnd("Tuple<Type, Type>");

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, Tuple<Type, Type>&)");

    return is;
}

template<class Type1, class Type2>
Ostream& operator<<(Ostream& os, const Tuple<Type1, Type2>& p)
{
    os  << token::BEGIN_LIST
        << p.first_ << token::SPACE
        << p.second_
        << token::END_LIST;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const Tuple<Type, Type>&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
