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
    Foam::functionObjectList

Description
    List of function objects with start(), execute() and end() functions
    that is called for each object.

See Also
    Foam::functionObject and Foam::OutputFilterFunctionObject

SourceFiles
    functionObjectList.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectList_H
#define functionObjectList_H

#include "PtrList.H"
#include "functionObject.H"
#include "SHA1Digest.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class functionObjectList Declaration
\*---------------------------------------------------------------------------*/

class functionObjectList
:
    private PtrList<functionObject>
{
    // Private data

        //- A list of SHA1 digests for the function object dictionaries
        List<SHA1Digest> digests_;

        //- Quick lookup of the index into functions/digests
        HashTable<label> indices_;

        const Time& time_;

        //- The parent dictionary containing a "functions" entry
        //  This entry can either be a list or a dictionary of
        //  functionObject specifications.
        const dictionary& parentDict_;

        //- Switch for the execution of the functionObjects
        bool execution_;

        //- Tracks if read() was called while execution is on
        bool updated_;


    // Private Member Functions

        //- Remove and return the function object pointer by name,
        //  and returns the old index via the parameter.
        //  Returns a NULL pointer (and index -1) if it didn't exist.
        functionObject* remove(const word&, label& oldIndex);

        //- Disallow default bitwise copy construct
        functionObjectList(const functionObjectList&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectList&);


public:

    // Constructors

        //- Construct from Time and the execution setting
        //  The functionObject specifications are read from the controlDict
        functionObjectList
        (
            const Time&,
            const bool execution=true
        );


        //- Construct from Time, a dictionary with "functions" entry
        //  and the execution setting.
        //  @param[in]  parentDict - the parent dictionary containing
        //    a "functions" entry, which can either be a list or a dictionary
        //    of functionObject specifications.
        functionObjectList
        (
            const Time&,
            const dictionary& parentDict,
            const bool execution=true
        );


    // Destructor

        virtual ~functionObjectList();


    // Member Functions

        //- Return the number of elements in the List.
        using PtrList<functionObject>::size;

        //- Return true if the List is empty (ie, size() is zero).
        using PtrList<functionObject>::empty;

        //- Access to the functionObjects
        using PtrList<functionObject>::operator[];

        //- Clear the list of function objects
        virtual void clear();


        //- Switch the function objects on
        virtual void on();

        //- Switch the function objects off
        virtual void off();

        //- Return the execution status (on/off) of the function objects
        virtual bool status() const;


        //- Called at the start of the time-loop
        virtual bool start();

        //- Called at each ++ or += of the time-loop
        virtual bool execute();

        //- Called when Time::run() determines that the time-loop exits
        virtual bool end();

        //- Read and set the function objects if their data have changed
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
