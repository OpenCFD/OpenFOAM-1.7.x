/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::validArgs.insert("file .. fileN");
    argList::validOptions.erase("case");
    argList::validOptions.insert("ext", "bak");

    argList args(argc, argv, false, true);

    if (args.additionalArgs().empty())
    {
        args.printUsage();
    }

    label ok = 0;

    forAll(args.additionalArgs(), argI)
    {
        const string& srcFile = args.additionalArgs()[argI];

        if (args.optionFound("ext"))
        {
            if (mvBak(srcFile, args.option("ext")))
            {
                ok++;
            }
        }
        else
        {
            if (mvBak(srcFile))
            {
                ok++;
            }
        }
    }

    Info<< "mvBak called for " << args.additionalArgs().size()
        << " files (moved " << ok << ")\n" << endl;

    return 0;
}


// ************************************************************************* //
