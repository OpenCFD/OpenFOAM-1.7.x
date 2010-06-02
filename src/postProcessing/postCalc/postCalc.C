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

Application
    postCalc

Description
    Generic wrapper for calculating a quantity at each time

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "timeSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    void tryCalc(const argList& args, const Time& runTime, const fvMesh& mesh)
    {
        FatalIOError.throwExceptions();

        try
        {
            calc(args, runTime, mesh);
        }
        catch(IOerror& err)
        {
            Warning<< err << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    Foam::argList::validOptions.insert("noWrite", "");
    Foam::argList::validOptions.insert("dict", "dictionary name");

#   include "setRootCase.H"
#   include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;

        mesh.readUpdate();

        Foam::tryCalc(args, runTime, mesh);

        Foam::Info<< Foam::endl;
    }

    return 0;
}


// ************************************************************************* //
