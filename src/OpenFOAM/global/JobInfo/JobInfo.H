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

Class
    Foam::JobInfo

Description
    Helper class for recording information about run/finished jobs.

    Writes the following files:
      - $FOAM_JOB_DIR/runningJobs
      - $FOAM_JOB_DIR/finishedJobs

SourceFiles
    JobInfo.C

\*---------------------------------------------------------------------------*/

#ifndef JobInfo_H
#define JobInfo_H

#include "dictionary.H"
#include "fileName.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JobInfo Declaration
\*---------------------------------------------------------------------------*/

class JobInfo
:
    public dictionary
{
    // Private data
        fileName runningJobPath_;
        fileName finishedJobPath_;
        cpuTime  cpuTime_;

    // Private member functions

        bool write(Ostream&) const;
        void end(const word& terminationType);


public:

    static bool constructed;
    static bool writeJobInfo;

    // Constructors

        //- Construct null
        JobInfo();


    // Destructor

        ~JobInfo();


    // Member Functions

        // Write

            void write() const;

            void end();

            void exit();

            void abort();

            void signalEnd() const;
};


extern JobInfo jobInfo;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
