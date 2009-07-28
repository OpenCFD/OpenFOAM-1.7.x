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
    Foam::BiIndirectList

Description
    Indexes into negList (negative index) or posList (zero or positive index).

SourceFiles
    BiIndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef BiIndirectList_H
#define BiIndirectList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BiIndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class BiIndirectList
{
    // Private data

        UList<T>& posList_;
        UList<T>& negList_;
        List<label> addressing_;


public:

    // Constructors

        //- Construct given the complete lists and the addressing array
        inline BiIndirectList
        (
            const UList<T>& posList,
            const UList<T>& negList,
            const UList<label>&
        );


    // Member Functions

        // Access

            inline label size() const;
            inline bool  empty() const;

            inline const UList<T>& posList() const;
            inline const UList<T>& negList() const;
            inline const List<label>& addressing() const;
            inline List<label>& addressing();

            //- Calculate index given whether index is into posList or negList
            inline static label posIndex(const label);
            inline static label negIndex(const label);


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BiIndirectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
