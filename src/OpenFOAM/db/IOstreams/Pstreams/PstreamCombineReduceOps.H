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

InClass
    Foam::Pstream

Description
    Combination-Reduction operation for a parallel run.  The
    information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes


\*---------------------------------------------------------------------------*/

#ifndef PstreamCombineReduceOps_H
#define PstreamCombineReduceOps_H

#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class T, class CombineOp>
void combineReduce
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    const CombineOp& cop
)
{
    Pstream::combineGather(comms, Value, cop);
    Pstream::combineScatter(comms, Value);
}


template <class T, class CombineOp>
void combineReduce(T& Value, const CombineOp& cop)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum)
    {
        Pstream::combineGather(Pstream::linearCommunication(), Value, cop);
        Pstream::combineScatter(Pstream::linearCommunication(), Value);
    }
    else
    {
        Pstream::combineGather(Pstream::treeCommunication(), Value, cop);
        Pstream::combineScatter(Pstream::treeCommunication(), Value);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
