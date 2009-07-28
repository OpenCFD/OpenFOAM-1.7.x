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

Description

\*---------------------------------------------------------------------------*/

#include "DynamicField.H"
#include "IOstreams.H"
#include "labelField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    {
        DynamicField<label> dl(10);
        Pout<< "null construct dl:" << dl << endl;
        dl.append(3);
        dl.append(2);
        dl.append(1);
        Pout<< "appending : dl:" << dl << endl;
    }

    {
        DynamicField<label> dl(IStringStream("(1 2 3)")());
        Pout<< "reading : dl:" << dl << endl;
    }

    {
        labelField lf(3);
        lf[0] = 1;
        lf[1] = 2;
        lf[2] = 3;
        DynamicField<label> dl;
        dl = lf;
        Pout<< "assigning from labelField : dl:" << dl << endl;
    }

    {
        labelField lf(3);
        lf[0] = 1;
        lf[1] = 2;
        lf[2] = 3;
        DynamicField<label> dl(lf);
        Pout<< "constructing from labelField dl:" << dl << endl;
    }


    Info<< "\nEnd\n";

    return 0;
}


// ************************************************************************* //
