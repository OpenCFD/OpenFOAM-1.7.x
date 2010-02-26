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

\*---------------------------------------------------------------------------*/

extern "C"
{
#include "scotch.h"
}

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static const char* notImplementedMessage =
"You are trying to use Scotch but do not have the Scotch library loaded.\n"
"This message is from the dummy Scotch stub library instead.\n"
"\n"
"Normally the Scotch library will be loaded through the LD_LIBRARY_PATH\n"
"environment variable but you are picking up this dummy library from the\n"
"$FOAM_LIBBIN/dummy directory. Please install Scotch and make sure the\n"
"libscotch.so is in your LD_LIBRARY_PATH.";


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_stratInit(SCOTCH_Strat *)
{
    Foam::FatalErrorIn("SCOTCH_stratInit(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_stratGraphMap(SCOTCH_Strat * const, const char *)
{
    Foam::FatalErrorIn("SCOTCH_stratGraphMap(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_stratSave
(
    const SCOTCH_Strat *,
    FILE * const
)
{
    Foam::FatalErrorIn("SCOTCH_stratSave(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
void SCOTCH_stratExit(SCOTCH_Strat *)
{
    Foam::FatalErrorIn("SCOTCH_stratExit(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}




#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphInit(SCOTCH_Graph *)
{
    Foam::FatalErrorIn("SCOTCH_graphInit(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphBuild
(
    SCOTCH_Graph*,
    const SCOTCH_Num,
    const SCOTCH_Num,
    const SCOTCH_Num*,
    const SCOTCH_Num*,
    const SCOTCH_Num*,
    const SCOTCH_Num*,
    const SCOTCH_Num,
    const SCOTCH_Num*,
    const SCOTCH_Num*
)
{
    Foam::FatalErrorIn("SCOTCH_graphBuild(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphCheck(const SCOTCH_Graph *)
{
    Foam::FatalErrorIn("SCOTCH_graphCheck(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
void SCOTCH_graphExit(SCOTCH_Graph *)
{
    Foam::FatalErrorIn("SCOTCH_graphExit(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}



#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_archInit(SCOTCH_Arch *)
{
    Foam::FatalErrorIn("SCOTCH_archInit(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_archCmpltw(SCOTCH_Arch*, const SCOTCH_Num, const SCOTCH_Num*)
{
    Foam::FatalErrorIn("SCOTCH_archCmpltw(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_archCmplt(SCOTCH_Arch *, const SCOTCH_Num)
{
    Foam::FatalErrorIn("SCOTCH_archCmplt(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
void SCOTCH_archExit(SCOTCH_Arch *)
{
    Foam::FatalErrorIn("SCOTCH_archExit(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}




#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphMapInit
(
    const SCOTCH_Graph*,
    SCOTCH_Mapping *,
    const SCOTCH_Arch *,
    SCOTCH_Num*
)
{
    Foam::FatalErrorIn("SCOTCH_graphMapInit(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphMapCompute
(
    const SCOTCH_Graph *,
    SCOTCH_Mapping *,
    const SCOTCH_Strat *
)
{
    Foam::FatalErrorIn("SCOTCH_graphMapCompute(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphMap
(
    const SCOTCH_Graph *,
    const SCOTCH_Arch *,
    const SCOTCH_Strat *,
    SCOTCH_Num *
)
{
    Foam::FatalErrorIn("SCOTCH_graphMap(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
    return -1;
}


#ifdef __cplusplus
extern "C"
#endif
void SCOTCH_graphMapExit (const SCOTCH_Graph *, SCOTCH_Mapping *)
{
    Foam::FatalErrorIn("SCOTCH_graphMapExit(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
