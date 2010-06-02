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
    foamUpgradeFvSolution

Description
    Simple tool to upgrade the syntax of system/fvSolution::solvers

Usage

    - foamUpgradeFvSolution [OPTION]

    @param -test \n
    Suppress writing the updated fvSolution file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "solution.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("test", "");

#   include "setRootCase.H"
#   include "createTime.H"

    IOdictionary solutionDict
    (
        IOobject
        (
            "fvSolution",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    label nChanged = 0;
    entry* e = solutionDict.lookupEntryPtr("solvers", false, false);
    if (e && e->isDict())
    {
        nChanged = solution::upgradeSolverDict(e->dict(), true);
    }

    Info<< nChanged << " solver settings changed" << nl << endl;
    if (nChanged)
    {
        if (args.optionFound("test"))
        {
            Info<< "-test option: no changes made" << nl << endl;
        }
        else
        {
            mv
            (
                solutionDict.objectPath(),
                solutionDict.objectPath() + ".old"
            );

            solutionDict.writeObject
            (
                IOstream::ASCII,
                IOstream::currentVersion,
                IOstream::UNCOMPRESSED
            );

            Info<< "Backup to    " << (solutionDict.objectPath() + ".old") << nl
                << "Write  to    " << solutionDict.objectPath() << nl << endl;
        }
    }

    return 0;
}


// ************************************************************************* //
