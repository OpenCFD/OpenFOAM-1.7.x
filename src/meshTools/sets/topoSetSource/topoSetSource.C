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

#include "topoSetSource.H"
#include "polyMesh.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(topoSetSource, 0);
defineRunTimeSelectionTable(topoSetSource, word);
defineRunTimeSelectionTable(topoSetSource, istream);

// Construct named object from dictionary
autoPtr<topoSetSource> topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_
            ->find(topoSetSourceType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "topoSetSource::New(const word&, "
            "const polyMesh&, const dictionary&)"
        )   << "Unknown topoSetSource type " << topoSetSourceType
            << endl << endl
            << "Valid topoSetSource types : " << endl
            << wordConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<topoSetSource>(cstrIter()(mesh, dict));
}


// Construct named object from Istream
autoPtr<topoSetSource> topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    Istream& is
)
{
    istreamConstructorTable::iterator cstrIter =
        istreamConstructorTablePtr_
            ->find(topoSetSourceType);

    if (cstrIter == istreamConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "topoSetSource::New(const word&, "
            "const polyMesh&, Istream&)"
        )   << "Unknown topoSetSource type " << topoSetSourceType
            << endl << endl
            << "Valid topoSetSource types : " << endl
            << istreamConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<topoSetSource>(cstrIter()(mesh, is));
}


} // End namespace Foam


Foam::HashTable<Foam::string>* Foam::topoSetSource::usageTablePtr_ = NULL;

template<>
const char* Foam::NamedEnum<Foam::topoSetSource::setAction, 8>::names[] =
{
    "clear",
    "new",
    "invert",
    "add",
    "delete",
    "subset",
    "list",
    "remove"
};


const Foam::NamedEnum<Foam::topoSetSource::setAction, 8>
    Foam::topoSetSource::actionNames_;


const Foam::string Foam::topoSetSource::illegalSource_
(
    "Illegal topoSetSource name"
);


Foam::Istream& Foam::topoSetSource::checkIs(Istream& is)
{
    if (is.good() && !is.eof())
    {
        return is;
    }
    else
    {
        FatalErrorIn("cellToFace::cellToFace") << "Istream not good"
            << exit(FatalError);

        return is;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::topoSetSource::addOrDelete
(
    topoSet& set,
    const label cellI,
    const bool add
) const
{
    if (add)
    {
        set.insert(cellI);
    }
    else
    {
        set.erase(cellI);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::topoSetSource::topoSetSource(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoSetSource::~topoSetSource()
{}


// ************************************************************************* //
