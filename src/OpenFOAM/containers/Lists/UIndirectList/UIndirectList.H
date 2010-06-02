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
    Foam::UIndirectList

Description
    A List with indirect addressing.
    Like IndirectList but does not store addressing.

SourceFiles
    UIndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef UIndirectList_H
#define UIndirectList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class T> class UIndirectList;
template<class T> Ostream& operator<<(Ostream&, const UIndirectList<T>&);

/*---------------------------------------------------------------------------*\
                        Class UIndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UIndirectList
{
    // Private data

        UList<T>& completeList_;
        const UList<label>& addressing_;


public:

    // Constructors

        //- Construct given the complete list and the addressing array
        inline UIndirectList(const UList<T>&, const UList<label>&);


    // Member Functions

        // Access

            inline label size() const;
            inline bool  empty() const;

            inline const UList<T>& completeList() const;
            inline const List<label>& addressing() const;


        // Member Operators

            //- Return the addressed elements as a List
            inline List<T> operator()() const;

            //- Return non-const access to an element
            inline T& operator[](const label);

            //- Return const access to an element
            inline const T& operator[](const label) const;

            //- Assignment from UList of addressed elements
            inline void operator=(const UList<T>&);

            //- Assignment of all entries to the given value
            inline void operator=(const T&);


    // Ostream operator

        //- Write UIndirectList to Ostream
        //  Binary output is currently still a bit of a problem
        friend Ostream& operator<< <T>
        (
            Ostream&,
            const UIndirectList<T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UIndirectListI.H"

#ifdef NoRepository
#   include "UIndirectListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
