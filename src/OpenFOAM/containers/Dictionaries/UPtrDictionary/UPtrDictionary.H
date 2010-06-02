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
    Foam::UPtrDictionary

Description
    Template dictionary class which does not manages the storage
    associated with it.

    It is derived from DictionaryBase instantiated on a non-memory managed
    form of intrusive doubly-linked list of \<T\>.

SourceFiles
    UPtrDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef UPtrDictionary_H
#define UPtrDictionary_H

#include "DictionaryBase.H"
#include "DLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class UPtrDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UPtrDictionary
:
    public DictionaryBase<DLList<T*>, T>
{

public:

    // Constructors

        //- Null constructor
        UPtrDictionary();

        //- Copy construct
        UPtrDictionary(const UPtrDictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "UPtrDictionary.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
