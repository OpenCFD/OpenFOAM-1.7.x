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

Class
    Foam::dlLibraryTable

Description
    A table of dynamically loaded libraries

SourceFiles
    dlLibraryTable.C

\*---------------------------------------------------------------------------*/

#ifndef dlLibraryTable_H
#define dlLibraryTable_H

#include "label.H"
#include "Hash.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dlLibraryTable Declaration
\*---------------------------------------------------------------------------*/

class dlLibraryTable
:
    public HashTable<fileName, void*, Hash<void*> >
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        dlLibraryTable(const dlLibraryTable&);

        //- Disallow default bitwise assignment
        void operator=(const dlLibraryTable&);


public:

    // Static data members

        //- Static data someStaticData
        static dlLibraryTable loadedLibraries;


    // Public classes

        //- Class whose construction causes the reading of dynamic libraries
        class readDlLibrary
        {
        public:

            //- Read all the libraries listed in the 'libsEntry' entry in the
            //  given dictionary if present
            readDlLibrary(const dictionary&, const word& libsEntry);
        };


    // Constructors

        //- Construct null
        dlLibraryTable();


    // Destructor

        ~dlLibraryTable();


    // Member Functions

        //- Open the named library
        static bool open(const fileName& name);

        //- Open all the libraries listed in the 'libsEntry' entry in the
        //  given dictionary if present
        static bool open(const dictionary&, const word& libsEntry);

        //- Open all the libraries listed in the 'libsEntry' entry in the
        //  given dictionary if present and check the additions
        //  to the give constructor table
        template<class TablePtr>
        static bool open
        (
            const dictionary&,
            const word& libsEntry,
            const TablePtr& tablePtr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "dlLibraryTableTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
