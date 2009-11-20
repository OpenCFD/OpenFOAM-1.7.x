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
    icoFoam

Description
    Incompressible laminar CFD code.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Perr<< "\nStarting transfers\n" << endl;

    vector data(0, 1, 2);

    if (Pstream::parRun())
    {
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            {
                Perr<< "slave sending to master "
                    << Pstream::masterNo() << endl;
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster << data;
            }

            Perr<< "slave receiving from master "
                << Pstream::masterNo() << endl;
            IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());
            fromMaster >> data;

            Perr<< data << endl;
        }
        else
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Perr << "master receiving from slave " << slave << endl;
                IPstream fromSlave(Pstream::scheduled, slave);
                fromSlave >> data;

                Perr<< data << endl;
            }

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Perr << "master sending to slave " << slave << endl;
                OPstream toSlave(Pstream::scheduled, slave);
                toSlave << data;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
