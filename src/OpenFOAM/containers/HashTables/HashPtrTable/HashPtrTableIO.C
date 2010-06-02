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

\*---------------------------------------------------------------------------*/

#include "HashPtrTable.H"
#include "Istream.H"
#include "Ostream.H"
#include "INew.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class INew>
void Foam::HashPtrTable<T, Key, Hash>::read(Istream& is, const INew& inewt)
{
    is.fatalCheck("HashPtrTable<T, Key, Hash>::read(Istream&, const INew&)");

    token firstToken(is);

    is.fatalCheck
    (
        "HashPtrTable<T, Key, Hash>::read(Istream&, const INew&) : "
        "reading first token"
    );

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        char delimiter = is.readBeginList("HashPtrTable<T, Key, Hash>");

        if (s)
        {
            if (2*s > this->tableSize_)
            {
                this->resize(2*s);
            }

            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<s; i++)
                {
                    Key key;
                    is >> key;
                    insert(key, inewt(key, is).ptr());

                    is.fatalCheck
                    (
                        "HashPtrTable<T, Key, Hash>::"
                        "read(Istream&, const INew&) : reading entry"
                    );
                }
            }
            else
            {
                FatalIOErrorIn
                (
                    "HashPtrTable<T, Key, Hash>::read(Istream&, const INew&)",
                    is
                )   << "incorrect first token, '(', found " << firstToken.info()
                    << exit(FatalIOError);
            }
        }

        // Read end of contents
        is.readEndList("HashPtrTable");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                "HashPtrTable<T, Key, Hash>::read(Istream&, const INew&)",
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);
            Key key;
            is >> key;
            insert(key, inewt(key, is).ptr());

            is.fatalCheck
            (
                "HashPtrTable<T, Key, Hash>::read(Istream&, const INew&) : "
                "reading entry"
            );

            is >> lastToken;
        }
    }
    else
    {
        FatalIOErrorIn
        (
            "HashPtrTable<T, Key, Hash>::read(Istream&, const INew&)",
            is
        )   << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    is.fatalCheck("HashPtrTable<T, Key, Hash>::read(Istream&, const INew&)");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class INew>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(Istream& is, const INew& inewt)
{
    read(is, inewt);
}


template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(Istream& is)
{
    read(is, INew<T>());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Istream& Foam::operator>>(Istream& is, HashPtrTable<T, Key, Hash>& L)
{
    L.clear();
    L.read(is, INew<T>());

    return is;
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const HashPtrTable<T, Key, Hash>& L
)
{
    // Write size and start delimiter
    os << nl << L.size() << nl << token::BEGIN_LIST << nl;

    // Write contents
    for
    (
        typename HashPtrTable<T, Key, Hash>::const_iterator iter = L.begin();
        iter != L.end();
        ++iter
    )
    {
        os << iter.key() << token::SPACE << *iter() << nl;
    }

    // Write end delimiter
    os << token::END_LIST;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const HashPtrTable&)");

    return os;
}


// ************************************************************************* //
