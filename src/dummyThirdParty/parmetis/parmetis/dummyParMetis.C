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
#include "parmetis.h"
}

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static const char* notImplementedMessage =
"You are trying to use Parmetis but do not have the Parmetis library loaded.\n"
"This message is from the dummy Parmetis stub library instead.\n"
"\n"
"Normally the Parmetis library will be loaded through the LD_LIBRARY_PATH\n"
"environment variable but you are picking up this dummy library from the\n"
"$FOAM_LIBBIN/dummy directory. Please install Parmetis and make sure the\n"
"libparmetis.so is in your LD_LIBRARY_PATH.";

#ifdef __cplusplus
extern "C"
#endif
void METIS_WPartGraphRecursive
(
    int *n,
    idxtype *xadj,
    idxtype *adjncy,
    idxtype *vwgt,
    idxtype *adjwgt,
    int *wgtflag,
    int *numflag,
    int *nparts,
    float *tpwgts,
    int *options,
    int *edgecut,
    idxtype *part
)
{
    Foam::FatalErrorIn("METIS_WPartGraphRecursive(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}


#ifdef __cplusplus
extern "C"
#endif
void METIS_PartGraphRecursive
(
    int *n,
    idxtype *xadj,
    idxtype *adjncy,
    idxtype *vwgt,
    idxtype *adjwgt,
    int *wgtflag,
    int *numflag,
    int *nparts,
    int *options,
    int *edgecut,
    idxtype *part
)
{
    Foam::FatalErrorIn("METIS_PartGraphRecursive(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}


#ifdef __cplusplus
extern "C"
#endif
void METIS_WPartGraphKway
(
    int *n,
    idxtype *xadj,
    idxtype *adjncy,
    idxtype *vwgt,
    idxtype *adjwgt,
    int *wgtflag,
    int *numflag,
    int *nparts,
    float *tpwgts,
    int *options,
    int *edgecut,
    idxtype *part
)
{
    Foam::FatalErrorIn("METIS_WPartGraphKway(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}


#ifdef __cplusplus
extern "C"
#endif
void METIS_PartGraphKway
(
    int *n,
    idxtype *xadj,
    idxtype *adjncy,
    idxtype *vwgt,
    idxtype *adjwgt,
    int *wgtflag,
    int *numflag,
    int *nparts,
    int *options,
    int *edgecut,
    idxtype *part
)
{
    Foam::FatalErrorIn("METIS_PartGraphKway(..)")
        << notImplementedMessage
        << Foam::exit(Foam::FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
