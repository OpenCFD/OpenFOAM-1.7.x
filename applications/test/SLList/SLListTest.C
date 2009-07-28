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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "IOstreams.H"
#include "SLList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    SLList<scalar> myList;

    for (int i = 0; i<10; i++)
    {
        myList.append(1.3*i);
    }

    myList.append(100.3);
    myList.append(500.3);

    Info<< nl << "And again using STL iterator: " << nl << endl;

    for
    (
        SLList<scalar>::iterator iter = myList.begin();
        iter != myList.end();
        ++iter
    )
    {
        Info<< "element:" << *iter << endl;
    }

    Info<< nl << "And again using STL const_iterator: " << nl << endl;

    const SLList<scalar>& const_myList = myList;

    for
    (
        SLList<scalar>::const_iterator iter2 = const_myList.begin();
        iter2 != const_myList.end();
        ++iter2
    )
    {
        Info<< "element:" << *iter2 << endl;
    }

    for
    (
        SLList<scalar>::iterator iter = myList.begin();
        iter != myList.end();
        ++iter
    )
    {
        Info<< "Removing element:" << *iter << endl;
        myList.remove(iter);
    }

    for
    (
        SLList<scalar>::const_iterator iter2 = const_myList.begin();
        iter2 != const_myList.end();
        ++iter2
    )
    {
        Info<< "element:" << *iter2 << endl;
    }



    for (int i = 0; i<10; i++)
    {
        myList.append(1.3*i);
    }

    myList.append(100.3);
    myList.append(500.3);

    Info<< nl << "Testing transfer: " << nl << endl;
    Info<< "original: " << myList << endl;

    SLList<scalar> newList;
    newList.transfer(myList);

    Info<< nl << "source: " << myList << nl
        << nl << "target: " << newList << endl;


    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
