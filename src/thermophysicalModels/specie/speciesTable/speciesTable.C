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

\*---------------------------------------------------------------------------*/

#include "speciesTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::speciesTable::setIndices()
{
    forAll (*this, i)
    {
        specieIndices_.insert(wordList::operator[](i), i);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from list of specie names
Foam::speciesTable::speciesTable(const wordList& specieNames)
:
    wordList(specieNames)
{
    setIndices();
}


// Construct from number of species and list of specie names
Foam::speciesTable::speciesTable(const label nSpecies, const char** specieNames)
:
    wordList(nSpecies)
{
    forAll (*this, i)
    {
        wordList::operator[](i) = specieNames[i];
    }

    setIndices();
}


// Construct from Istream
Foam::speciesTable::speciesTable(Istream& is)
:
    wordList(is)
{
    setIndices();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, speciesTable& st)
{
    is >> static_cast<wordList&>(st);
    st.setIndices();

    return is;
}

// ************************************************************************* //
