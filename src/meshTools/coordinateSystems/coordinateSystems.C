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

#include "coordinateSystems.H"
#include "IOPtrList.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystems, 0);
    defineTemplateTypeNameAndDebug(IOPtrList<coordinateSystem>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::coordinateSystems(const IOobject& io)
:
    IOPtrList<coordinateSystem>(io)
{}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    const PtrList<coordinateSystem>& lst
)
:
    IOPtrList<coordinateSystem>(io, lst)
{}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    const Xfer<PtrList<coordinateSystem> >& lst
)
:
    IOPtrList<coordinateSystem>(io, lst)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// Read construct from registry, or return previously registered
const Foam::coordinateSystems& Foam::coordinateSystems::New
(
    const objectRegistry& obr
)
{
    if (obr.foundObject<coordinateSystems>(typeName))
    {
        return obr.lookupObject<coordinateSystems>(typeName);
    }
    else
    {
        return obr.store
        (
            new coordinateSystems
            (
                IOobject
                (
                    typeName,
                    "constant",
                    obr,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::coordinateSystems::find(const word& keyword) const
{
    forAll(*this, i)
    {
        if (keyword == operator[](i).name())
        {
            return i;
        }
    }

    return -1;
}


bool Foam::coordinateSystems::found(const word& keyword) const
{
    return find(keyword) >= 0;
}


Foam::wordList Foam::coordinateSystems::toc() const
{
    wordList keywords(size());

    forAll(*this, i)
    {
        keywords[i] = operator[](i).name();
    }

    return keywords;
}


bool Foam::coordinateSystems::writeData(Ostream& os, bool subDict) const
{
    os << nl << size() << nl << token::BEGIN_LIST;

    forAll(*this, i)
    {
        os << nl;
        operator[](i).writeDict(os, subDict);
    }

    os << token::END_LIST << nl;

    return os.good();
}


// ************************************************************************* //
