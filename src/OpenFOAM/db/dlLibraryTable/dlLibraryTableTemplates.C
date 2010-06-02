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
#include "dictionary.H"
#include "fileNameList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TablePtr>
bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry,
    const TablePtr& tablePtr
)
{
    if (dict.found(libsEntry))
    {
        fileNameList libNames(dict.lookup(libsEntry));

        bool allOpened = (libNames.size() > 0);

        forAll(libNames, i)
        {
            const fileName& libName = libNames[i];

            label nEntries = 0;

            if (tablePtr)
            {
                nEntries = tablePtr->size();
            }

            bool opened = dlLibraryTable::open(libName);
            allOpened = opened && allOpened;

            if (opened && (!tablePtr || tablePtr->size() <= nEntries))
            {
                WarningIn
                (
                    "dlLibraryTable::open"
                    "(const dictionary& dict, const word& libsEntry, "
                    "const TablePtr tablePtr)"
                )   << "library " << libName
                    << " did not introduce any new entries"
                    << endl << endl;
            }
        }

        return allOpened;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
