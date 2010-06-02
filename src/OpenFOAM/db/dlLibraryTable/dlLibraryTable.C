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

#include "dlLibraryTable.H"

#include <dlfcn.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::dlLibraryTable Foam::dlLibraryTable::loadedLibraries;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dlLibraryTable::dlLibraryTable()
:
    HashTable<fileName, void*, Hash<void*> >()
{}


Foam::dlLibraryTable::readDlLibrary::readDlLibrary
(
    const dictionary& dict,
    const word& libsEntry
)
{
    open(dict, libsEntry);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dlLibraryTable::~dlLibraryTable()
{
    forAllIter(dlLibraryTable, *this, iter)
    {
        dlclose(iter.key());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dlLibraryTable::open(const fileName& functionLibName)
{
    if (functionLibName.size())
    {
        void* functionLibPtr = 
            dlopen(functionLibName.c_str(), RTLD_LAZY|RTLD_GLOBAL);

        if (!functionLibPtr)
        {
            WarningIn
            (
                "dlLibraryTable::open(const fileName& functionLibName)"
            )   << "could not load " << dlerror()
                << endl;

            return false;
        }
        else
        {
            if (!loadedLibraries.found(functionLibPtr))
            {
                loadedLibraries.insert(functionLibPtr, functionLibName);
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        return false;
    }
}


bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry
)
{
    if (dict.found(libsEntry))
    {
        fileNameList libNames(dict.lookup(libsEntry));

        bool allOpened = (libNames.size() > 0);

        forAll(libNames, i)
        {
            allOpened = dlLibraryTable::open(libNames[i]) && allOpened;
        }

        return allOpened;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
