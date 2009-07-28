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
    Foam::Map

Description
    A HashTable to objects of type \<T\> with a label key.

See Also
    PtrMap

\*---------------------------------------------------------------------------*/

#ifndef Map_H
#define Map_H

#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Map Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Map
:
    public HashTable<T, label, Hash<label> >
{

public:

    typedef typename HashTable<T, label, Hash<label> >::iterator iterator;

    typedef typename HashTable<T, label, Hash<label> >::const_iterator
        const_iterator;

    // Constructors

        //- Construct given initial size
        Map(const label size = 128)
        :
            HashTable<T, label, Hash<label> >(size)
        {}

        //- Construct from Istream
        Map(Istream& is)
        :
            HashTable<T, label, Hash<label> >(is)
        {}

        //- Construct as copy
        Map(const Map<T>& map)
        :
            HashTable<T, label, Hash<label> >(map)
        {}

        //- Construct by transferring the parameter contents
        Map(const Xfer<Map<T> >& map)
        :
            HashTable<T, label, Hash<label> >(map)
        {}

        //- Construct by transferring the parameter contents
        Map(const Xfer<HashTable<T, label, Hash<label> > >& map)
        :
            HashTable<T, label, Hash<label> >(map)
        {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
