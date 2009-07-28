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
    foamCalc

Description
    Generic wrapper for calculating a quantity at each time. Split into four
    phases:
        1. Intialise
        2. Pre-time calculation loop
        3. Calculation loop
        4. Post-calculation loop

\*---------------------------------------------------------------------------*/

#include "timeSelector.H"
#include "calcType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    Foam::argList::validOptions.insert("noWrite", "");
    Foam::argList::validOptions.insert("dict", "dictionary name");

    if (argc < 2)
    {
        FatalError
            << "No utility has been supplied" << nl
            << exit(FatalError);
    }

    word utilityName = argv[1];

    Foam::autoPtr<Foam::calcType> utility
    (
        calcType::New(utilityName)
    );

    utility().tryInit();

#   include "setRootCase.H"
#   include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
#   include "createMesh.H"

    utility().tryPreCalc(args, runTime, mesh);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;

        mesh.readUpdate();

        utility().tryCalc(args, runTime, mesh);

        Foam::Info<< Foam::endl;
    }

    utility().tryPostCalc(args, runTime, mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
