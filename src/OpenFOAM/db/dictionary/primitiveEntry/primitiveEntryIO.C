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

Description
    PrimitiveEntry constructor from Istream and Ostream output operator.

\*---------------------------------------------------------------------------*/

#include "primitiveEntry.H"
#include "OSspecific.H"
#include "functionEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::primitiveEntry::append
(
    const token& currToken,
    const dictionary& dict,
    Istream& is
)
{
    if (currToken.isWord())
    {
        const word& w = currToken.wordToken();

        if
        (
            w.size() == 1
         || (
                !(w[0] == '$' && expandVariable(w, dict))
             && !(w[0] == '#' && expandFunction(w, dict, is))
            )
        )
        {
            newElmt(tokenIndex()++) = currToken;
        }
    }
    else
    {
        newElmt(tokenIndex()++) = currToken;
    }
}


void Foam::primitiveEntry::append(const tokenList& varTokens)
{
    forAll(varTokens, i)
    {
        newElmt(tokenIndex()++) = varTokens[i];
    }
}


bool Foam::primitiveEntry::expandVariable
(
    const word& w,
    const dictionary& dict
)
{
    word varName = w(1, w.size()-1);

    // lookup the variable name in the given dictionary....
    const entry* ePtr = dict.lookupEntryPtr(varName, true, true);

    // ...if defined insert its tokens into this
    if (ePtr != NULL)
    {
        append(ePtr->stream());
        return true;
    }
    else
    {
        // if not in the dictionary see if it is an environment
        // variable

        string enVarString = getEnv(varName);

        if (enVarString.size())
        {
            append(tokenList(IStringStream('(' + enVarString + ')')()));
            return true;
        }

        return false;
    }
}


bool Foam::primitiveEntry::expandFunction
(
    const word& keyword,
    const dictionary& parentDict,
    Istream& is
)
{
    word functionName = keyword(1, keyword.size()-1);
    return functionEntry::execute(functionName, parentDict, *this, is);
}


bool Foam::primitiveEntry::read(const dictionary& dict, Istream& is)
{
    is.fatalCheck
    (
        "primitiveEntry::readData(const dictionary&, Istream&)"
    );

    label blockCount = 0;
    token currToken;

    if
    (
        !is.read(currToken).bad()
     && currToken.good()
     && currToken != token::END_STATEMENT
    )
    {
        append(currToken, dict, is);

        if
        (
            currToken == token::BEGIN_BLOCK
         || currToken == token::BEGIN_LIST
        )
        {
            blockCount++;
        }

        while
        (
            !is.read(currToken).bad()
         && currToken.good()
         && !(currToken == token::END_STATEMENT && blockCount == 0)
        )
        {
            if
            (
                currToken == token::BEGIN_BLOCK
             || currToken == token::BEGIN_LIST
            )
            {
                blockCount++;
            }
            else if
            (
                currToken == token::END_BLOCK
             || currToken == token::END_LIST
            )
            {
                blockCount--;
            }

            append(currToken, dict, is);
        }
    }

    is.fatalCheck
    (
        "primitiveEntry::readData(const dictionary&, Istream&)"
    );

    if (currToken.good())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::primitiveEntry::readEntry(const dictionary& dict, Istream& is)
{
    label keywordLineNumber = is.lineNumber();
    tokenIndex() = 0;

    if (read(dict, is))
    {
        setSize(tokenIndex());
        tokenIndex() = 0;
    }
    else
    {
        FatalIOErrorIn
        (
            "primitiveEntry::readEntry(const dictionary&, Istream&)",
            is
        )   << "ill defined primitiveEntry starting at keyword '"
            << keyword() << '\''
            << " on line " << keywordLineNumber
            << " and ending at line " << is.lineNumber()
            << exit(FatalIOError);
    }
}


Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const dictionary& dict,
    Istream& is
)
:
    entry(key),
    ITstream
    (
        is.name() + "::" + key,
        tokenList(10),
        is.format(),
        is.version()
    )
{
    readEntry(dict, is);
}


Foam::primitiveEntry::primitiveEntry(const keyType& key, Istream& is)
:
    entry(key),
    ITstream
    (
        is.name() + "::" + key,
        tokenList(10),
        is.format(),
        is.version()
    )
{
    readEntry(dictionary::null, is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::primitiveEntry::write(Ostream& os) const
{
    os.writeKeyword(keyword());

    for (label i=0; i<size(); i++)
    {
        os << operator[](i);

        if (i < size()-1)
        {
            os << token::SPACE;
        }
    }

    os << token::END_STATEMENT << endl;
}


// * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<primitiveEntry>& ip
)
{
    const primitiveEntry& e = ip.t_;

    e.print(os);

    const label nPrintTokens = 10;

    os  << "    primitiveEntry '" << e.keyword() << "' comprises ";

    for (label i=0; i<min(e.size(), nPrintTokens); i++)
    {
        os  << nl << "        " << e[i].info();
    }

    if (e.size() > nPrintTokens)
    {
        os  << " ...";
    }

    os  << endl;

    return os;
}


// ************************************************************************* //
