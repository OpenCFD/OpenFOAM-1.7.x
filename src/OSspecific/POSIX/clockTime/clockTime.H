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
    Foam::clockTime

Description
    Starts timing (using rtc) and returns elapsed time from start.
    Better resolution (2uSec instead of ~20mSec) than cpuTime.

SourceFiles
    clockTime.C

\*---------------------------------------------------------------------------*/

#ifndef clockTime_H
#define clockTime_H

#include <sys/types.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clockTime Declaration
\*---------------------------------------------------------------------------*/

class clockTime
{
    // Private data

        struct timeval startTime_;
        mutable struct timeval lastTime_;
        mutable struct timeval newTime_;

        static void getTime(struct timeval& t);

        static double timeDifference
        (
            const struct timeval& start,
            const struct timeval& end
        );


public:

    // Constructors

        //- Construct from components
        clockTime();


    // Member Functions

        // Access

            //- Returns CPU time from start of run
            double elapsedTime() const;

            //- Returns CPU time from last call of clockTimeIncrement()
            double timeIncrement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
