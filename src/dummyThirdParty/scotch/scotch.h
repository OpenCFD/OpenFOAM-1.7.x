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

Namespace
    scotch

Description
    Dummy stub for scotch library functions. Only implements the absolute
    minimum we are using.

SourceFiles
    scotch.C

\*---------------------------------------------------------------------------*/

#ifndef scotch_H
#define scotch_H

#ifdef __cplusplus
extern "C"
#endif
{
#include <stdio.h>
}


#ifndef SCOTCH_Strat
#define SCOTCH_Strat int
#define SCOTCH_Num int
#define SCOTCH_Graph int
#define SCOTCH_Arch int
#define SCOTCH_Mapping int
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class scotch Declaration
\*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_stratInit(SCOTCH_Strat *);

#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_stratGraphMap(SCOTCH_Strat * const, const char *);

#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_stratSave
(
    const SCOTCH_Strat *,
    FILE * const
);

#ifdef __cplusplus
extern "C"
#endif
void SCOTCH_stratExit(SCOTCH_Strat *);



#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphInit(SCOTCH_Graph *);

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
);

#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphCheck(const SCOTCH_Graph *);

#ifdef __cplusplus
extern "C"
#endif
void SCOTCH_graphExit(SCOTCH_Graph *);



#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_archInit(SCOTCH_Arch *);

#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_archCmpltw(SCOTCH_Arch*, const SCOTCH_Num, const SCOTCH_Num* const);

#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_archCmplt(SCOTCH_Arch *, const SCOTCH_Num);

#ifdef __cplusplus
extern "C"
#endif
void SCOTCH_archExit(SCOTCH_Arch *);



#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphMapInit
(
    const SCOTCH_Graph*,
    SCOTCH_Mapping *,
    const SCOTCH_Arch *,
    SCOTCH_Num*
);

#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphMapCompute
(
    const SCOTCH_Graph *,
    SCOTCH_Mapping *,
    const SCOTCH_Strat *
);

#ifdef __cplusplus
extern "C"
#endif
int SCOTCH_graphMap
(
    const SCOTCH_Graph *,
    const SCOTCH_Arch *,
    const SCOTCH_Strat *,
    SCOTCH_Num *
);

#ifdef __cplusplus
extern "C"
#endif
void SCOTCH_graphMapExit (const SCOTCH_Graph *, SCOTCH_Mapping *);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
