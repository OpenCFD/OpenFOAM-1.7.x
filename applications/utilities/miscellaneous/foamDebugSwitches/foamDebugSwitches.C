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
    Write out all library debug switches

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "IFstream.H"
#include "IOobject.H"
#include "HashSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("new", "");
    argList::validOptions.insert("old", "");

    Foam::argList args(argc, argv);

    wordList currDebug(debug::debugSwitches().toc());
    wordList currInfo(debug::infoSwitches().toc());
    wordList currOpt(debug::optimisationSwitches().toc());

    if (args.optionFound("old") || args.optionFound("new"))
    {
        dictionary controlDict(IFstream(findEtcFile("controlDict", true))());

        wordHashSet oldDebug
        (
            controlDict.subDict("DebugSwitches").toc()
        );

        wordHashSet oldInfo
        (
            controlDict.subDict("InfoSwitches").toc()
        );

        wordHashSet oldOpt
        (
            controlDict.subDict("OptimisationSwitches").toc()
        );


        wordHashSet hashset;
        wordList listing;


        // list old switches - but this can't work since the (old) inserted
        // switches are in both sets
        // Workaround:
        //  1. run without any options (get complete list)
        //  2. comment out DebugSwitches, run again with -new to find new ones
        //     and do a diff
        if (args.optionFound("old"))
        {
            IOobject::writeDivider(Info);

            hashset = oldDebug;
            hashset -= currDebug;
            listing = hashset.toc();
            sort(listing);
            Info<< "old DebugSwitches: " << listing << endl;

            hashset = oldInfo;
            hashset -= currInfo;
            listing = hashset.toc();
            sort(listing);
            Info<< "old InfoSwitches: " << listing << endl;

            hashset = oldOpt;
            hashset -= currOpt;
            listing = hashset.toc();
            sort(listing);
            Info<< "old OptimisationSwitches: " << listing << endl;
        }

        // list new switches
        if (args.optionFound("new"))
        {
            IOobject::writeDivider(Info);

            hashset = currDebug;
            hashset -= oldDebug;

            listing = hashset.toc();
            sort(listing);
            Info<< "new DebugSwitches: " << listing << endl;

            hashset = currInfo;
            hashset -= oldInfo;
            listing = hashset.toc();
            sort(listing);
            Info<< "new InfoSwitches: " << listing << endl;

            hashset = currOpt;
            hashset -= oldOpt;
            listing = hashset.toc();
            sort(listing);
            Info<< "new OptimisationSwitches: " << listing << endl;
        }
    }
    else
    {
        IOobject::writeDivider(Info);

        sort(currDebug);
        Info<< "DebugSwitches: " << currDebug << endl;

        sort(currInfo);
        Info<< "InfoSwitches: " << currInfo << endl;

        sort(currOpt);
        Info<< "OptimisationSwitches: " << currOpt << endl;
    }



    Info<< "done" << endl;

    return 0;
}


// ************************************************************************* //
