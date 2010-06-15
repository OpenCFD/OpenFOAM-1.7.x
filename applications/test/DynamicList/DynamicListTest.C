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

#include "DynamicList.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<DynamicList<label, 1, 0> > ldl(2);

    ldl[0](0) = 0;
    ldl[0](2) = 2;
    ldl[0](3) = 3;
    ldl[0](1) = 1;

    ldl[0].setCapacity(5);    // increase allocated size
    ldl[1].setCapacity(10);   // increase allocated size
    ldl[0].reserve(15);       // should increase allocated size
    ldl[1].reserve(5);        // should not decrease allocated size
    ldl[1](3) = 2;            // allocates space and sets value

    // this works without a segfault, but doesn't change the list size
    ldl[0][4] = 4;

    ldl[1] = 3;

    Info<< "<ldl>" << ldl << "</ldl>" << nl << "sizes: ";
    forAll(ldl, i)
    {
        Info<< " " << ldl[i].size() << "/" << ldl[i].capacity();
    }
    Info<< endl;

    List<List<label> > ll(2);
    ll[0].transfer(ldl[0]);
    ll[1].transfer(ldl[1].shrink());

    Info<< "<ldl>" << ldl << "</ldl>" << nl << "sizes: ";
    forAll(ldl, i)
    {
        Info<< " " << ldl[i].size() << "/" << ldl[i].capacity();
    }
    Info<< endl;

    Info<< "<ll>" << ll << "</ll>" << nl << endl;


    // test the transfer between DynamicLists
    DynamicList<label, 1, 0> dlA;
    DynamicList<label, 1, 0> dlB;

    for (label i = 0; i < 5; i++)
    {
        dlA.append(i);
    }
    dlA.setCapacity(10);

    Info<< "<dlA>" << dlA << "</dlA>" << nl << "sizes: "
        << " " << dlA.size() << "/" << dlA.capacity() << endl;

    dlB.transfer(dlA);

    // provokes memory error if previous transfer did not maintain
    // the correct allocated space
    dlB[6] = 6;

    Info<< "Transferred to dlB" << endl;
    Info<< "<dlA>" << dlA << "</dlA>" << nl << "sizes: "
        << " " << dlA.size() << "/" << dlA.capacity() << endl;
    Info<< "<dlB>" << dlB << "</dlB>" << nl << "sizes: "
        << " " << dlB.size() << "/" << dlB.capacity() << endl;

    // try with a normal list:
    List<label> lstA;
    lstA.transfer(dlB);
    Info<< "Transferred to normal list" << endl;
    Info<< "<lstA>" << lstA << "</lstA>" << nl << "sizes: "
        << " " << lstA.size() << endl;
    Info<< "<dlB>" << dlB << "</dlB>" << nl << "sizes: "
        << " " << dlB.size() << "/" << dlB.capacity() << endl;

    // Copy back and append a few time
    for (label i=0; i < 3; i++)
    {
        dlB.append(lstA);
    }

    Info<< "appended list a few times" << endl;
    Info<< "<dlB>" << dlB << "</dlB>" << nl << "sizes: "
        << " " << dlB.size() << "/" << dlB.capacity() << endl;

    // assign the list (should maintain allocated space)
    dlB = lstA;
    Info<< "assigned list" << endl;
    Info<< "<dlB>" << dlB << "</dlB>" << nl << "sizes: "
        << " " << dlB.size() << "/" << dlB.capacity() << endl;


    // Copy back and append a few time
    for (label i=0; i < 3; i++)
    {
        dlB.append(lstA);
    }


    // check allocation granularity
    DynamicList<label, 6, 0> dlC;

    Info<< "<dlC>" << dlC << "</dlC>" << nl << "sizes: "
        << " " << dlC.size() << "/" << dlC.capacity() << endl;

    dlC.reserve(dlB.size());
    dlC = dlB;

    Info<< "<dlC>" << dlC << "</dlC>" << nl << "sizes: "
        << " " << dlC.size() << "/" << dlC.capacity() << endl;

    List<label> lstB(dlC.xfer());

    Info<< "Transferred to normal list via the xfer() method" << endl;
    Info<< "<lstB>" << lstB << "</lstB>" << nl << "sizes: "
        << " " << lstB.size() << endl;
    Info<< "<dlC>" << dlC << "</dlC>" << nl << "sizes: "
        << " " << dlC.size() << "/" << dlC.capacity() << endl;

    DynamicList<label> dlD(lstB.xfer());

    Info<< "Transfer construct from normal list" << endl;
    Info<< "<lstB>" << lstB << "</lstB>" << nl << "sizes: "
        << " " << lstB.size() << endl;
    Info<< "<dlD>" << dlD << "</dlD>" << nl << "sizes: "
        << " " << dlD.size() << "/" << dlD.capacity() << endl;

    DynamicList<label,10> dlE1(10);
    DynamicList<label> dlE2(dlE1);   // construct dissimilar

    Info<< "<dlE1>" << dlE1 << "</dlE1>" << nl << "sizes: "
        << " " << dlE1.size() << "/" << dlE1.capacity() << endl;
    Info<< "<dlE2>" << dlE2 << "</dlE2>" << nl << "sizes: "
        << " " << dlE2.size() << "/" << dlE2.capacity() << endl;

    for (label elemI=0; elemI < 5; ++elemI)
    {
        dlE1.append(4 - elemI);
        dlE2.append(elemI);
    }

    Info<< "<dlE2>" << dlE2 << "</dlE2>" << endl;

    DynamicList<label> dlE3(dlE2);   // construct identical
    Info<< "<dlE3>" << dlE3 << "</dlE3>" << endl;

    dlE3 = dlE1;   // assign dissimilar
    Info<< "<dlE3>" << dlE3 << "</dlE3>" << endl;

    dlE3 = dlE2;   // assign identical
    Info<< "<dlE3>" << dlE3 << "</dlE3>" << endl;


    Info<< "\nEnd\n";

    return 0;
}


// ************************************************************************* //
