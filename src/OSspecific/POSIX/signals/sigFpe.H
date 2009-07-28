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
    Foam::sigFpe

Description
    Set up trapping for floating point exceptions (signal FPE).

    Controlled by two env vars:
    @param FOAM_SIGFPE \n
        exception trapping
    @param FOAM_SETNAN \n
        initialization of all malloced memory to NaN. If FOAM_SIGFPE
        also set, this will cause usage of uninitialized scalars to trigger
        an abort.

SourceFiles
    sigFpe.C

\*---------------------------------------------------------------------------*/

#ifndef sigFpe_H
#define sigFpe_H

#include "OSspecific.H"
#include <signal.h>

#if defined(linux) || defined(linuxAMD64) || defined(linuxIA64)
#    define LINUX
#endif

#if defined(LINUX) && defined(__GNUC__)
#    define LINUX_GNUC
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sigFpe Declaration
\*---------------------------------------------------------------------------*/

class sigFpe
{
    // Private data

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;

#       ifdef LINUX

            //- Saved old malloc
            static void *(*old_malloc_hook)(size_t, const void *);

            //- nan malloc function. From malloc_hook manpage.
            static void* my_malloc_hook(size_t size, const void *caller);

#       endif


    // Static data members

#   ifdef LINUX_GNUC

        //- Handler for caught signals
        static void sigFpeHandler(int);

#   endif


public:


    // Constructors

        sigFpe();


    // Destructor

        ~sigFpe();


    // Member functions

        void set(const bool verbose);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
