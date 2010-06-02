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

Description
    Interrogates a case and prints information to screen

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dictionary.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("times", "");
    argList::validOptions.insert("dictionary", "dictionary name");
    argList::validOptions.insert("keywords", "");
    argList::validOptions.insert("entry", "entry name");

#   include "setRootCase.H"

    Info<< endl;

    if (args.optionFound("times"))
    {
        instantList times
        (
            Foam::Time::findTimes(args.rootPath()/args.caseName())
        );

        forAll (times, i)
        {
            Info<< times[i].name() << endl;
        }
    }

    if (args.optionFound("dictionary"))
    {
        fileName dictFileName
        (
            args.rootPath()/args.caseName()/args.option("dictionary")
        );

        IFstream dictFile(dictFileName);

        if (dictFile.good())
        {
            dictionary dict(dictFile);

            if (args.optionFound("keywords") && !args.optionFound("entry"))
            {
                for
                (
                    IDLList<entry>::iterator iter = dict.begin();
                    iter != dict.end();
                    ++iter
                )
                {
                    Info<< iter().keyword() << endl;
                }
            }
            else if (args.optionFound("entry"))
            {
                wordList entryNames
                (
                    fileName(args.option("entry")).components(':')
                );

                if (dict.found(entryNames[0]))
                {
                    const entry* entPtr = &dict.lookupEntry
                    (
                        entryNames[0],
                        false,
                        true            // wildcards
                    );

                    for (int i=1; i<entryNames.size(); i++)
                    {
                        if (entPtr->dict().found(entryNames[i]))
                        {
                            entPtr = &entPtr->dict().lookupEntry
                            (
                                entryNames[i],
                                false,
                                true    // wildcards
                            );
                        }
                        else
                        {
                            FatalErrorIn(args.executable())
                                << "Cannot find sub-entry " << entryNames[i]
                                << " in entry " << args.option("entry")
                                << " in dictionary " << dictFileName;
                            FatalError.exit(3);
                        }
                    }

                    if (args.optionFound("keywords"))
                    {
                        /*
                        if (ent[1] != token::BEGIN_BLOCK)
                        {
                            FatalErrorIn(args.executable())
                                << "Cannot find entry "
                                << args.option("entry")
                                << " in dictionary " << dictFileName
                                << " is not a sub-dictionary";
                            FatalError.exit(4);
                        }
                        */

                        const dictionary& dict(entPtr->dict());
                        for
                        (
                            IDLList<entry>::const_iterator iter = dict.begin();
                            iter != dict.end();
                            ++iter
                        )
                        {
                            Info<< iter().keyword() << endl;
                        }
                    }
                    else
                    {
                        Info<< *entPtr << endl;
                    }
                }
                else
                {
                    FatalErrorIn(args.executable())
                        << "Cannot find entry "
                        << entryNames[0]
                        << " in dictionary " << dictFileName;
                    FatalError.exit(2);
                }
            }
            else
            {
                Info<< dict;
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Cannot open file " << dictFileName;
            FatalError.exit(1);
        }
    }

    return 0;
}


// ************************************************************************* //
