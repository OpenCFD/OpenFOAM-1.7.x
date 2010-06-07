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

#include "IndirectList.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<double> completeList(10);

    forAll(completeList, i)
    {
        completeList[i] = 0.1*i;
    }

    List<label> addresses(5);
    addresses[0] = 1;
    addresses[1] = 0;
    addresses[2] = 7;
    addresses[3] = 8;
    addresses[4] = 5;

    IndirectList<double> idl(completeList, addresses);

    forAll(idl, i)
    {
        Info<< idl[i] << token::SPACE;
    }

    Info<< endl;

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
