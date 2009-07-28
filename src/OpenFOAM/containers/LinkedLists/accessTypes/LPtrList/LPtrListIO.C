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

\*---------------------------------------------------------------------------*/

#include "LPtrList.H"
#include "Istream.H"
#include "Ostream.H"
#include "INew.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class LListBase, class T>
template<class INew>
void Foam::LPtrList<LListBase, T>::read(Istream& is, const INew& iNew)
{
    is.fatalCheck
    (
        "LPtrList<LListBase, T>::read(Istream&, const INew&)"
    );

    token firstToken(is);

    is.fatalCheck
    (
        "LPtrList<LListBase, T>::read(Istream&, const INew&) : "
        "reading first token"
    );

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        char delimiter = is.readBeginList("LPtrList<LListBase, T>");

        if (s)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<s; i++)
                {
                    append(iNew(is).ptr());

                    is.fatalCheck
                    (
                        "LPtrList<LListBase, T>::read(Istream&, const INew&) : "
                        "reading entry"
                    );
                }
            }
            else
            {
                T* tPtr = iNew(is).ptr();
                append(tPtr);

                is.fatalCheck
                (
                    "LPtrList<LListBase, T>::read(Istream&, const INew&) : "
                    "reading entry"
                );

                for (label i=1; i<s; i++)
                {
                    append(tPtr->clone().ptr());
                }
            }
        }

        // Read end of contents
        is.readEndList("LPtrList<LListBase, T>");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                "LPtrList<LListBase, T>::read(Istream&, const INew&)",
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        is.fatalCheck("LPtrList<LListBase, T>::read(Istream&, const INew&)");

        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);
            append(iNew(is).ptr());

            is >> lastToken;
            is.fatalCheck
            (
                "LPtrList<LListBase, T>::read(Istream&, const INew&)"
            );
        }
    }
    else
    {
        FatalIOErrorIn
        (
            "LPtrList<LListBase, T>::read(Istream&, const INew&)",
            is
        )   << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    is.fatalCheck("LPtrList<LListBase, T>::read(Istream&, const INew&)");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
template<class INew>
Foam::LPtrList<LListBase, T>::LPtrList(Istream& is, const INew& iNew)
{
    read(is, iNew);
}


template<class LListBase, class T>
Foam::LPtrList<LListBase, T>::LPtrList(Istream& is)
{
    read(is, INew<T>());
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Istream& Foam::operator>>(Istream& is, LPtrList<LListBase, T>& L)
{
    L.clear();
    L.read(is, INew<T>());

    return is;
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const LPtrList<LListBase, T>& lst)
{
    // Write size
    os << nl << lst.size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST << nl;

    // Write contents
    for
    (
        typename LPtrList<LListBase, T>::const_iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        os << iter() << nl;
    }

    // Write end of contents
    os << token::END_LIST;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const LPtrList<LListBase, T>&)");

    return os;
}

// ************************************************************************* //
