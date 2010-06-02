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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "wordReList.H"

#include "IOstreams.H"
#include "IStringStream.H"
#include "scalar.H"
#include "vector.H"
#include "ListOps.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("reList", "reList");
    argList::validOptions.insert("wordList", "wordList");
    argList::validOptions.insert("stringList", "stringList");
    argList::validOptions.insert("float", "xx");
    argList::validOptions.insert("flag", "");

#   include "setRootCase.H"

    List<vector> list1(IStringStream("1 ((0 1 2))")());
    Info<< "list1: " << list1 << endl;

    List<vector> list2(IStringStream("((0 1 2) (3 4 5) (6 7 8))")());
    Info<< "list2: " << list2 << endl;

    list1.append(list2);
    Info<< "list1.append(list2): " << list1 << endl;

    Info<< findIndex(list2, vector(3, 4, 5)) << endl;

    list2.setSize(10, vector(1, 2, 3));
    Info<< "list2: " << list2 << endl;

    List<vector> list3(list2.xfer());
    Info<< "Transferred via the xfer() method" << endl;
    Info<< "list2: " << list2 << nl
        << "list3: " << list3 << endl;


    // Subset
    const labelList map(IStringStream("2 (0 2)")());
    List<vector> subList3(list3, map);
    Info<< "Elements " << map << " out of " << list3
        << " => " << subList3 << endl;

    wordReList reLst;
    wordList wLst;
    stringList sLst;


    scalar xxx(-1);

    if (args.optionFound("flag"))
    {
        Info<<"-flag:" << args.option("flag") << endl;
    }

    if (args.optionReadIfPresent<scalar>("float", xxx))
    {
        Info<<"read float " << xxx << endl;
    }

    if (args.optionFound("reList"))
    {
        reLst = args.optionReadList<wordRe>("reList");
    }

    if (args.optionFound("wordList"))
    {
        wLst = args.optionReadList<word>("wordList");
    }

    if (args.optionFound("stringList"))
    {
        sLst = args.optionReadList<string>("stringList");
    }

    Info<< nl
        << "-reList: " << reLst << nl
        << "-wordList: " << wLst << nl
        << "-stringList: " << sLst << endl;

    return 0;
}

// ************************************************************************* //
