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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a single value constructed from the values
    on individual processors using a user-specified operator.

\*---------------------------------------------------------------------------*/

#include "OPstream.H"
#include "IPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T, class BinaryOp>
void Pstream::gather
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    const BinaryOp& bop
)
{
    if (Pstream::parRun())
    {
        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo()];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            T value;

            if (contiguous<T>())
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    myComm.below()[belowI],
                    reinterpret_cast<char*>(&value),
                    sizeof(T)
                );
            }
            else
            {
                IPstream fromBelow(Pstream::scheduled, myComm.below()[belowI]);
                fromBelow >> value;
            }

            Value = bop(Value, value);
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            if (contiguous<T>())
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T)
                );
            }
            else
            {
                OPstream toAbove(Pstream::scheduled, myComm.above());
                toAbove << Value;
            }
        }
    }
}


template <class T, class BinaryOp>
void Pstream::gather(T& Value, const BinaryOp& bop)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum)
    {
        gather(Pstream::linearCommunication(), Value, bop);
    }
    else
    {
        gather(Pstream::treeCommunication(), Value, bop);
    }
}


template <class T>
void Pstream::scatter(const List<Pstream::commsStruct>& comms, T& Value)
{
    if (Pstream::parRun())
    {
        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            if (contiguous<T>())
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    myComm.above(),
                    reinterpret_cast<char*>(&Value),
                    sizeof(T)
                );
            }
            else
            {
                IPstream fromAbove(Pstream::scheduled, myComm.above());
                fromAbove >> Value;
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            if (contiguous<T>())
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    myComm.below()[belowI],
                    reinterpret_cast<const char*>(&Value),
                    sizeof(T)
                );
            }
            else
            {
                OPstream toBelow(Pstream::scheduled,myComm.below()[belowI]);
                toBelow << Value;
            }
        }
    }
}


template <class T>
void Pstream::scatter(T& Value)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum)
    {
        scatter(Pstream::linearCommunication(), Value);
    }
    else
    {
        scatter(Pstream::treeCommunication(), Value);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
