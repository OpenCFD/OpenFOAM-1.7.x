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
    surfaceMeshImport

Description
    Import from various third-party surface formats into surfMesh
    with optional scaling or transformations (rotate/translate)
    on a coordinateSystem.

Usage
    - surfaceMeshImport inputFile [OPTION]

    @param -clean \n
    Perform some surface checking/cleanup on the input surface.

    @param -name \<name\> \n
    Specify an alternative surface name when writing.

    @param -scaleIn \<scale\> \n
    Specify a scaling factor when reading files.

    @param -scaleOut \<scale\> \n
    Specify a scaling factor when writing files.

    @param -dict \<dictionary\> \n
    Specify an alternative dictionary for constant/coordinateSystems.

    @param -from \<coordinateSystem\> \n
    Specify a coordinate System when reading files.

    @param -to \<coordinateSystem\> \n
    Specify a coordinate System when writing files.

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"

#include "MeshedSurfaces.H"
#include "coordinateSystems.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("inputFile");
    argList::validOptions.insert("name",  "name");
    argList::validOptions.insert("clean", "");
    argList::validOptions.insert("scaleIn",  "scale");
    argList::validOptions.insert("scaleOut", "scale");
    argList::validOptions.insert("dict", "coordinateSystemsDict");
    argList::validOptions.insert("from", "sourceCoordinateSystem");
    argList::validOptions.insert("to",   "targetCoordinateSystem");

#   include "setRootCase.H"
#   include "createTime.H"

    const stringList& params = args.additionalArgs();

    // try for the latestTime, but create "constant" as needed
    instantList Times = runTime.times();
    if (Times.size())
    {
        label startTime = Times.size()-1;
        runTime.setTime(Times[startTime], startTime);
    }
    else
    {
        runTime.setTime(instant(0, runTime.constant()), 0);
    }


    fileName importName(params[0]);
    word exportName("default");
    args.optionReadIfPresent("name", exportName);

    // check that reading is supported
    if (!MeshedSurface<face>::canRead(importName, true))
    {
        return 1;
    }


    // get the coordinate transformations
    autoPtr<coordinateSystem> fromCsys;
    autoPtr<coordinateSystem> toCsys;

    if (args.optionFound("from") || args.optionFound("to"))
    {
        autoPtr<IOobject> ioPtr;

        if (args.optionFound("dict"))
        {
            fileName dictPath(args.option("dict"));

            ioPtr.set
            (
                new IOobject
                (
                    (
                        isDir(dictPath)
                      ? dictPath/coordinateSystems::typeName
                      : dictPath
                    ),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
        }
        else
        {
            ioPtr.set
            (
                new IOobject
                (
                    coordinateSystems::typeName,
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
        }


        if (!ioPtr->headerOk())
        {
            FatalErrorIn(args.executable())
                << "Cannot open coordinateSystems file\n    "
                << ioPtr->objectPath() << nl
                << exit(FatalError);
        }

        coordinateSystems csLst(ioPtr());

        if (args.optionFound("from"))
        {
            const word csName(args.option("from"));

            label csId = csLst.find(csName);
            if (csId < 0)
            {
                FatalErrorIn(args.executable())
                    << "Cannot find -from " << csName << nl
                    << "available coordinateSystems: " << csLst.toc() << nl
                    << exit(FatalError);
            }

            fromCsys.reset(new coordinateSystem(csLst[csId]));
        }

        if (args.optionFound("to"))
        {
            const word csName(args.option("to"));

            label csId = csLst.find(csName);
            if (csId < 0)
            {
                FatalErrorIn(args.executable())
                    << "Cannot find -to " << csName << nl
                    << "available coordinateSystems: " << csLst.toc() << nl
                    << exit(FatalError);
            }

            toCsys.reset(new coordinateSystem(csLst[csId]));
        }


        // maybe fix this later
        if (fromCsys.valid() && toCsys.valid())
        {
            FatalErrorIn(args.executable())
                << "Only allowed  '-from' or '-to' option at the moment."
                << exit(FatalError);
        }
    }



    MeshedSurface<face> surf(importName);

    if (args.optionFound("clean"))
    {
        surf.cleanup(true);
    }


    scalar scaleIn = 0;
    if (args.optionReadIfPresent("scaleIn", scaleIn) && scaleIn > 0)
    {
        Info<< " -scaleIn " << scaleIn << endl;
        surf.scalePoints(scaleIn);
    }

    if (fromCsys.valid())
    {
        Info<< " -from " << fromCsys().name() << endl;
        tmp<pointField> tpf = fromCsys().localPosition(surf.points());
        surf.movePoints(tpf());
    }

    if (toCsys.valid())
    {
        Info<< " -to " << toCsys().name() << endl;
        tmp<pointField> tpf = toCsys().globalPosition(surf.points());
        surf.movePoints(tpf());
    }

    scalar scaleOut = 0;
    if (args.optionReadIfPresent("scaleOut", scaleOut) && scaleOut > 0)
    {
        Info<< " -scaleOut " << scaleOut << endl;
        surf.scalePoints(scaleOut);
    }

    surfMesh smesh
    (
        IOobject
        (
            exportName,
            runTime.constant(),
            runTime
        ),
        surf.xfer()
    );


    Info<< "writing surfMesh:\n  " << smesh.objectPath() << endl;
    smesh.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
