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

Class
    Foam::ParRunControl

Description
    Helper class for initializing parallel jobs from the command arguments.

\*---------------------------------------------------------------------------*/

#ifndef parRun_H
#define parRun_H

#include "OPstream.H"
#include "IPstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ParRunControl Declaration
\*---------------------------------------------------------------------------*/

class ParRunControl
{
    bool RunPar;

public:

    ParRunControl()
    :
        RunPar(false)
    {}

    ~ParRunControl()
    {
        if (RunPar)
        {
            Info<< "Finalising parallel run" << endl;
            Pstream::exit(0);
        }
    }

    void runPar(int& argc, char**& argv)
    {
        RunPar = true;

        if (!Pstream::init(argc, argv))
        {
            Info<< "Failed to start parallel run" << endl;
            Pstream::exit(1);
        }
    }

    bool parRun() const
    {
        return RunPar;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
