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

\*---------------------------------------------------------------------------*/

#include "NamedEnum.H"
#include "stringList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Enum, int nEnum>
Foam::NamedEnum<Enum, nEnum>::NamedEnum()
:
    HashTable<int>(2*nEnum)
{
    for (int i=0; i<nEnum; i++)
    {
        if (!names[i] || names[i][0] == '\0')
        {
            stringList goodNames(i);

            for (int j = 0; j < i; j++)
            {
                goodNames[j] = names[j];
            }

            FatalErrorIn("NamedEnum<Enum, nEnum>::NamedEnum()")
                << "Illegal enumeration name at position " << i << endl
                << "after entries " << goodNames << ".\n"
                << "Possibly your NamedEnum<Enum, nEnum>::names array"
                << " is not of size " << nEnum << endl
                << abort(FatalError);
        }
        insert(names[i], i);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Enum, int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::read(Istream& is) const
{
    word name(is);

    HashTable<int>::const_iterator iter = find(name);

    if (iter == HashTable<int>::end())
    {
        FatalIOErrorIn
        (
            "NamedEnum<Enum, nEnum>::read(Istream&) const", is
        )   << name << " is not in enumeration: "
            << toc() << exit(FatalIOError);
    }

    return Enum(iter());
}


template<class Enum, int nEnum>
void Foam::NamedEnum<Enum, nEnum>::write(const Enum e, Ostream& os) const
{
    os  << operator[](e);
}


// ************************************************************************* //
