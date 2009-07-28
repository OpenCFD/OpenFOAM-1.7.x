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

Application
    surfaceMeshExport

Description
    Export from surfMesh to various third-party surface formats with
    optional scaling or transformations (rotate/translate) on a
    coordinateSystem.

Usage
    - surfaceMeshExport outputFile [OPTION]

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
    argList::validArgs.append("outputFile");
    argList::validOptions.insert("name",  "name");
    argList::validOptions.insert("clean", "");
    argList::validOptions.insert("scaleIn",  "scale");
    argList::validOptions.insert("scaleOut", "scale");
    argList::validOptions.insert("dict", "coordinateSystemsDict");
    argList::validOptions.insert("from", "sourceCoordinateSystem");
    argList::validOptions.insert("to",   "targetCoordinateSystem");

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());
    const stringList& params = args.additionalArgs();

    fileName exportName(params[0]);
    word importName("default");
    args.optionReadIfPresent("name", importName);

    // check that writing is supported
    if (!MeshedSurface<face>::canWriteType(exportName.ext(), true))
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


    surfMesh smesh
    (
        IOobject
        (
            importName,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "read surfMesh:\n  " << smesh.objectPath() << endl;


    // Simply copy for now, but really should have a separate write method

    MeshedSurface<face> surf(smesh);

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


    surf.writeStats(Info);
    Info<< endl;

    Info<< "writing " << exportName << endl;
    surf.write(exportName);

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
