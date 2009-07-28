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

#include "mpi.h"

#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"

#include <cstring>
#include <cstdlib>
#include <csignal>

#if defined(WM_SP)
#   define MPI_SCALAR MPI_FLOAT
#elif defined(WM_DP)
#   define MPI_SCALAR MPI_DOUBLE
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::Pstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("GAMMANP", "number of instances");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::Pstream::init(int& argc, char**& argv)
{
    MPI_Init(&argc, &argv);

    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcNo_);

    if (numprocs <= 1)
    {
        FatalErrorIn("Pstream::init(int& argc, char**& argv)")
            << "bool Pstream::init(int& argc, char**& argv) : "
               "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }

    procIDs_.setSize(numprocs);

    forAll(procIDs_, procNo)
    {
        procIDs_[procNo] = procNo;
    }

    setParRun();

#   ifndef SGIMPI
    string bufferSizeName = getEnv("MPI_BUFFER_SIZE");

    if (bufferSizeName.size())
    {
        int bufferSize = atoi(bufferSizeName.c_str());

        if (bufferSize)
        {
            MPI_Buffer_attach(new char[bufferSize], bufferSize);
        }
    }
    else
    {
        FatalErrorIn("Pstream::init(int& argc, char**& argv)")
            << "Pstream::init(int& argc, char**& argv) : "
            << "environment variable MPI_BUFFER_SIZE not defined"
            << Foam::abort(FatalError);
    }
#   endif

    int processorNameLen;
    char processorName[MPI_MAX_PROCESSOR_NAME];

    MPI_Get_processor_name(processorName, &processorNameLen);

    //signal(SIGABRT, stop);

    // Now that nprocs is known construct communication tables.
    initCommunicationSchedule();

    return true;
}


void Foam::Pstream::exit(int errnum)
{
#   ifndef SGIMPI
    int size;
    char* buff;
    MPI_Buffer_detach(&buff, &size);
    delete[] buff;
#   endif

    if (errnum == 0)
    {
        MPI_Finalize();
        ::exit(errnum);
    }
    else
    {
        MPI_Abort(MPI_COMM_WORLD, errnum);
    }
}


void Foam::Pstream::abort()
{
    MPI_Abort(MPI_COMM_WORLD, 1);
}


void Foam::reduce(scalar& Value, const sumOp<scalar>& bop)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (Pstream::nProcs() <= Pstream::nProcsSimpleSum)
    {
        if (Pstream::master())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                scalar value;

                if
                (
                    MPI_Recv
                    (
                        &value,
                        1,
                        MPI_SCALAR,
                        Pstream::procID(slave),
                        Pstream::msgType(),
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    )
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "MPI_Recv failed"
                        << Foam::abort(FatalError);
                }

                Value = bop(Value, value);
            }
        }
        else
        {
            if
            (
                MPI_Send
                (
                    &Value,
                    1,
                    MPI_SCALAR,
                    Pstream::procID(Pstream::masterNo()),
                    Pstream::msgType(),
                    MPI_COMM_WORLD
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "MPI_Send failed"
                    << Foam::abort(FatalError);
            }
        }


        if (Pstream::master())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                if
                (
                    MPI_Send
                    (
                        &Value,
                        1,
                        MPI_SCALAR,
                        Pstream::procID(slave),
                        Pstream::msgType(),
                        MPI_COMM_WORLD
                    )
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "MPI_Send failed"
                        << Foam::abort(FatalError);
                }
            }
        }
        else
        {
            if
            (
                MPI_Recv
                (
                    &Value,
                    1,
                    MPI_SCALAR,
                    Pstream::procID(Pstream::masterNo()),
                    Pstream::msgType(),
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "MPI_Recv failed"
                    << Foam::abort(FatalError);
            }
        }
    }
    else
    {
        scalar sum;
        MPI_Allreduce(&Value, &sum, 1, MPI_SCALAR, MPI_SUM, MPI_COMM_WORLD);
        Value = sum;

        /*
        int myProcNo = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();

        //
        // receive from children
        //
        int level = 1;
        int thisLevelOffset = 2;
        int childLevelOffset = thisLevelOffset/2;
        int childProcId = 0;

        while
        (
            (childLevelOffset < nProcs)
         && (myProcNo % thisLevelOffset) == 0
        )
        {
            childProcId = myProcNo + childLevelOffset;

            scalar value;

            if (childProcId < nProcs)
            {
                if
                (
                    MPI_Recv
                    (
                        &value,
                        1,
                        MPI_SCALAR,
                        Pstream::procID(childProcId),
                        Pstream::msgType(),
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    )
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "MPI_Recv failed"
                        << Foam::abort(FatalError);
                }

	            Value = bop(Value, value);
	        }

            level++;
            thisLevelOffset <<= 1;
            childLevelOffset = thisLevelOffset/2;
        }

        //
        // send and receive from parent
        //
        if (!Pstream::master())
        {
            int parentId = myProcNo - (myProcNo % thisLevelOffset);

            if
            (
                MPI_Send
                (
                    &Value,
                    1,
                    MPI_SCALAR,
                    Pstream::procID(parentId),
                    Pstream::msgType(),
                    MPI_COMM_WORLD
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "MPI_Send failed"
                    << Foam::abort(FatalError);
            }

            if
            (
                MPI_Recv
                (
                    &Value,
                    1,
                    MPI_SCALAR,
                    Pstream::procID(parentId),
                    Pstream::msgType(),
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                )
            )
            {
                FatalErrorIn
                (
                    "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                )   << "MPI_Recv failed"
                    << Foam::abort(FatalError);
            }
        }


        //
        // distribute to my children
        //
        level--;
        thisLevelOffset >>= 1;
        childLevelOffset = thisLevelOffset/2;

        while (level > 0)
        {
            childProcId = myProcNo + childLevelOffset;

            if (childProcId < nProcs)
            {
                if
                (
                    MPI_Send
                    (
                        &Value,
                        1,
                        MPI_SCALAR,
                        Pstream::procID(childProcId),
                        Pstream::msgType(),
                        MPI_COMM_WORLD
                    )
                )
                {
                    FatalErrorIn
                    (
                        "reduce(scalar& Value, const sumOp<scalar>& sumOp)"
                    )   << "MPI_Send failed"
                        << Foam::abort(FatalError);
                }
            }

            level--;
            thisLevelOffset >>= 1;
            childLevelOffset = thisLevelOffset/2;
        }
        */
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
