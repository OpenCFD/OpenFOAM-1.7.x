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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "sigFpe.H"

#include "JobInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"

#ifdef LINUX_GNUC

#   ifndef __USE_GNU
#       define __USE_GNU
#   endif

#   include <fenv.h>
#   include <malloc.h>

#elif defined(sgiN32) || defined(sgiN32Gcc)

#   include <sigfpe.h>

#endif

#include <stdint.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

struct sigaction Foam::sigFpe::oldAction_;


#if defined(LINUX)

void *(*Foam::sigFpe::old_malloc_hook)(size_t, const void *) = NULL;

void* Foam::sigFpe::my_malloc_hook(size_t size, const void *caller)
{
    void *result;

    // Restore all old hooks
    __malloc_hook = old_malloc_hook;

    // Call recursively
    result = malloc (size);

    // initialize to signalling nan
#   ifdef WM_SP

    const uint32_t sNAN = 0x7ff7fffflu;

    int nScalars = size / sizeof(scalar);

    uint32_t* dPtr = reinterpret_cast<uint32_t*>(result);

    for (int i = 0; i < nScalars; i++)
    {
        *dPtr++ = sNAN;
    }

#   else

    const uint64_t sNAN = 0x7ff7ffffffffffffllu;

    int nScalars = size/sizeof(scalar);

    uint64_t* dPtr = reinterpret_cast<uint64_t*>(result);

    for (int i = 0; i < nScalars; i++)
    {
        *dPtr++ = sNAN;
    }

#   endif

    // Restore our own hooks
    __malloc_hook = my_malloc_hook;

    return result;
}

#endif


#ifdef LINUX_GNUC

void Foam::sigFpe::sigFpeHandler(int)
{
    // Reset old handling
    if (sigaction(SIGFPE, &oldAction_, NULL) < 0)
    {
        FatalErrorIn
        (
            "Foam::sigSegv::sigFpeHandler()"
        )   << "Cannot reset SIGFPE trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo.signalEnd();

    error::printStack(Perr);

    // Throw signal (to old handler)
    raise(SIGFPE);
}

#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    oldAction_.sa_handler = NULL;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    if (env("FOAM_SIGFPE"))
    {
#       ifdef LINUX_GNUC

        // Reset signal
        if (oldAction_.sa_handler && sigaction(SIGFPE, &oldAction_, NULL) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigFpe::~sigFpe()"
            )   << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }

#       endif
    }

    if (env("FOAM_SETNAN"))
    {
#       ifdef LINUX_GNUC

        // Reset to standard malloc
        if (oldAction_.sa_handler)
        {
            __malloc_hook = old_malloc_hook;
        }

#       endif
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigFpe::set(const bool verbose)
{
    if (oldAction_.sa_handler)
    {
        FatalErrorIn
        (
            "Foam::sigFpe::set()"
        )   << "Cannot call sigFpe::set() more than once"
            << abort(FatalError);
    }

    if (env("FOAM_SIGFPE"))
    {
        if (verbose)
        {
            Info<< "SigFpe : Enabling floating point exception trapping"
                << " (FOAM_SIGFPE)." << endl;
        }

#       ifdef LINUX_GNUC

        feenableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        struct sigaction newAction;
        newAction.sa_handler = sigFpeHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(SIGFPE, &newAction, &oldAction_) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigFpe::set()"
            )   << "Cannot set SIGFPE trapping"
                << abort(FatalError);
        }


#       elif defined(sgiN32) || defined(sgiN32Gcc)

        sigfpe_[_DIVZERO].abort=1;
        sigfpe_[_OVERFL].abort=1;
        sigfpe_[_INVALID].abort=1;

        sigfpe_[_DIVZERO].trace=1;
        sigfpe_[_OVERFL].trace=1;
        sigfpe_[_INVALID].trace=1;

        handle_sigfpes
        (
            _ON,
            _EN_DIVZERO
          | _EN_INVALID
          | _EN_OVERFL,
            0,
            _ABORT_ON_ERROR,
            NULL
        );

#       endif
    }


    if (env("FOAM_SETNAN"))
    {
        if (verbose)
        {
            Info<< "SetNaN : Initialising allocated memory to NaN"
                << " (FOAM_SETNAN)." << endl;
        }

#       ifdef LINUX_GNUC

        // Set our malloc
        __malloc_hook = Foam::sigFpe::my_malloc_hook;

#       endif
    }
}


// ************************************************************************* //
