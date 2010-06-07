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

Application
    surfaceConvert

Description
    Converts from one surface mesh format to another

Usage
    - surfaceConvert inputFile outputFile [OPTION]

    @param -clean \n
    Perform some surface checking/cleanup on the input surface

    @param -scale \<scale\> \n
    Specify a scaling factor for writing the files

    @param -group \n
    Orders faces by region

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "OSspecific.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("inputFile");
    argList::validArgs.append("outputFile");
    argList::validOptions.insert("clean", "");
    argList::validOptions.insert("scale", "scale");
    argList::validOptions.insert("group", "");

    argList args(argc, argv);
    const stringList& params = args.additionalArgs();

    fileName importName(params[0]);
    fileName exportName(params[1]);

    if (importName == exportName)
    {
        FatalErrorIn(args.executable())
            << "Output file " << exportName << " would overwrite input file."
            << exit(FatalError);
    }

    Info<< "Reading : " << importName << endl;
    triSurface surf(importName);

    Info<< "Read surface:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    if (args.optionFound("clean"))
    {
        Info<< "Cleaning up surface" << endl;
        surf.cleanup(true);

        Info<< "After cleaning up surface:" << endl;
        surf.writeStats(Info);
        Info<< endl;
    }

    bool sortByRegion = args.optionFound("group");

    if (sortByRegion)
    {
        Info<< "Reordering faces into groups; one per region." << endl;
    }
    else
    {
        Info<< "Maintaining face ordering" << endl;
    }

    Info<< "writing " << exportName;

    scalar scaleFactor = 0;
    if (args.optionReadIfPresent("scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< " with scaling " << scaleFactor;
        surf.scalePoints(scaleFactor);
    }
    Info<< endl;

    surf.write(exportName, sortByRegion);

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
