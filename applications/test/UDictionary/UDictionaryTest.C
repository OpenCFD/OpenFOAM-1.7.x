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

Application

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "IOstreams.H"
#include "UDictionary.H"

using namespace Foam;

class ent
:
    public UDictionary<ent>::link
{
    word keyword_;
    int i_;

public:

    ent(const word& keyword, int i)
    :
        keyword_(keyword),
        i_(i)
    {}

    const word& keyword() const
    {
        return keyword_;
    }

    friend Ostream& operator<<(Ostream& os, const ent& e)
    {
        os << e.keyword_ << ' ' << e.i_ << endl;
        return os;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    UDictionary<ent>* dictPtr = new UDictionary<ent>;
    UDictionary<ent>& dict = *dictPtr;

    for (int i = 0; i<10; i++)
    {
        ent* ePtr = new ent(word("ent") + name(i), i);
        dict.append(ePtr->keyword(), ePtr);
        dict.swapUp(ePtr);
    }

    Info<< dict << endl;

    dict.swapDown(dict.first());

    for
    (
        UDictionary<ent>::const_iterator iter = dict.begin();
        iter != dict.end();
        ++iter
    )
    {
        Info<< "element : " << *iter;
    }

    Info<< dict.toc() << endl;

    delete dictPtr;

    dictPtr = new UDictionary<ent>;
    UDictionary<ent>& dict2 = *dictPtr;

    for (int i = 0; i<10; i++)
    {
        ent* ePtr = new ent(word("ent") + name(i), i);
        dict2.append(ePtr->keyword(), ePtr);
        dict2.swapUp(ePtr);
    }

    Info<< dict2 << endl;

    dict2.remove("ent9");
    dict2.UILList<DLListBase, ent>::remove(dict2.first());

    Info<< dict2 << endl;


    Info<< nl << "Testing transfer: " << nl << endl;
    Info<< "original: " << dict2 << endl;

    UDictionary<ent> newDict;
    newDict.transfer(dict2);

    Info<< nl << "source: " << dict2 << nl
        << "keys: " << dict2.toc() << nl
        << "target: " << newDict << nl
        << "keys: " << newDict.toc() << endl;

    Info<< nl << "Done." << endl;

    return 0;
}


// ************************************************************************* //
