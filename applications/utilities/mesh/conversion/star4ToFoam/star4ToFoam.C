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
    star4ToFoam

Description
    Converts a Star-CD (v4) pro-STAR mesh into OpenFOAM format.

Usage
    - star4ToFoam [OPTION] ccmMesh\n
      convert pro-STAR mesh to OpenFOAM

    @param -ascii \n
    Write in ASCII format instead of binary

    @param -scale \<factor\>\n
    Specify an alternative geometry scaling factor.
    The default is @b 0.001 (scale @em [mm] to @em [m]).

    @param -solids \n
    Treat any solid cells present just like fluid cells.
    The default is to discard them.

Note
    - baffles are written as interfaces for later use

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "STARCDMeshReader.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("pro-STAR prefix");
    argList::validOptions.insert("ascii", "");
    argList::validOptions.insert("scale", "scale");
    argList::validOptions.insert("solids", "");

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());
    const stringList& params = args.additionalArgs();

    // default rescale from [mm] to [m]
    scalar scaleFactor = 0.001;
    if (args.optionReadIfPresent("scale", scaleFactor))
    {
        if (scaleFactor <= 0)
        {
            scaleFactor = 1;
        }
    }

    if (args.optionFound("solids"))
    {
        meshReaders::STARCD::keepSolids = true;
    }

    // default to binary output, unless otherwise specified
    IOstream::streamFormat format = IOstream::BINARY;
    if (args.optionFound("ascii"))
    {
        format = IOstream::ASCII;
    }

    // increase the precision of the points data
    IOstream::defaultPrecision(10);

    // remove extensions and/or trailing '.'
    fileName prefix = fileName(params[0]).lessExt();

    meshReaders::STARCD reader(prefix, runTime, scaleFactor);

    autoPtr<polyMesh> mesh = reader.mesh(runTime);
    reader.writeMesh(mesh, format);


    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
