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
    C linkage

Description
    Dummy stub for metis library functions. From Metis-4.0 documentation.
    Only implements the absolute minimum we are using.

SourceFiles
    dummyMetis.C

\*---------------------------------------------------------------------------*/

#ifndef metis_H
#define metis_H

#ifndef idxtype
#define idxtype int
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class metis Declaration
\*---------------------------------------------------------------------------*/

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
);

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
);

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
);

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
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
