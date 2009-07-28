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

#include "error.H"
#include "OStringStream.H"
#include "fileName.H"
#include "dictionary.H"
#include "JobInfo.H"
#include "Pstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::error::error(const string& title)
:
    std::exception(),
    messageStream(title, messageStream::FATAL),
    functionName_("unknown"),
    sourceFileName_("unknown"),
    sourceFileLineNumber_(0),
    abort_(env("FOAM_ABORT")),
    throwExceptions_(false),
    messageStreamPtr_(new OStringStream())
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::error(const string& title) : cannot open error stream"
            << endl;
        exit(1);
    }
}


Foam::error::error(const dictionary& errDict)
:
    std::exception(),
    messageStream(errDict),
    functionName_(errDict.lookup("functionName")),
    sourceFileName_(errDict.lookup("sourceFileName")),
    sourceFileLineNumber_(readLabel(errDict.lookup("sourceFileLineNumber"))),
    abort_(env("FOAM_ABORT")),
    throwExceptions_(false),
    messageStreamPtr_(new OStringStream())
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::error(const dictionary& errDict) : "
               "cannot open error stream"
            << endl;
        exit(1);
    }
}


Foam::error::error(const error& err)
:
    std::exception(),
    messageStream(err),
    functionName_(err.functionName_),
    sourceFileName_(err.sourceFileName_),
    sourceFileLineNumber_(err.sourceFileLineNumber_),
    abort_(err.abort_),
    throwExceptions_(err.throwExceptions_),
    messageStreamPtr_(new OStringStream(*err.messageStreamPtr_))
{
    //*messageStreamPtr_ << err.message();
}


Foam::error::~error() throw()
{
    delete messageStreamPtr_;
}


Foam::OSstream& Foam::error::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    messageStreamPtr_->rewind();
    functionName_ = functionName;
    sourceFileName_ = sourceFileName;
    sourceFileLineNumber_ = sourceFileLineNumber;

    return operator OSstream&();
}


Foam::OSstream& Foam::error::operator()
(
    const string& functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    return operator()
    (
        functionName.c_str(),
        sourceFileName,
        sourceFileLineNumber
    );
}


Foam::error::operator OSstream&()
{
    if (!messageStreamPtr_->good())
    {
        Perr<< endl
            << "error::operator OSstream&() : error stream has failed"
            << endl;
        printStack(Perr);
        abort();
    }

    return *messageStreamPtr_;
}


Foam::error::operator dictionary() const
{
    dictionary errDict;

    string oneLineMessage(message());
    oneLineMessage.replaceAll('\n', ' ');

    errDict.add("type", word("Foam::error"));
    errDict.add("message", oneLineMessage);
    errDict.add("function", functionName());
    errDict.add("sourceFile", sourceFileName());
    errDict.add("sourceFileLineNumber", sourceFileLineNumber());

    return errDict;
}


Foam::string Foam::error::message() const
{
    return messageStreamPtr_->str();
}


void Foam::error::exit(const int errNo)
{
    if (!throwExceptions_ && JobInfo::constructed)
    {
        jobInfo.add("FatalError", operator dictionary());
        jobInfo.exit();
    }

    if (abort_)
    {
        printStack(*this);
        Perr<< endl << *this << endl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
        abort();
    }

    if (Pstream::parRun())
    {
        Perr<< endl << *this << endl
            << "\nFOAM parallel run exiting\n" << endl;
        Pstream::exit(errNo);
    }
    else
    {
        if (throwExceptions_)
        {
            throw *this;
        }
        else
        {
            Perr<< endl << *this << endl
                << "\nFOAM exiting\n" << endl;
            ::exit(1);
        }
    }
}


void Foam::error::abort()
{
    if (!throwExceptions_ && JobInfo::constructed)
    {
        jobInfo.add("FatalError", operator dictionary());
        jobInfo.abort();
    }

    if (abort_)
    {
        printStack(*this);
        Perr<< endl << *this << endl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
        ::abort();
    }

    if (Pstream::parRun())
    {
        printStack(*this);
        Perr<< endl << *this << endl
            << "\nFOAM parallel run aborting\n" << endl;
        Pstream::abort();
    }
    else
    {
        if (throwExceptions_)
        {
            throw *this;
        }
        else
        {
            printStack(*this);
            Perr<< endl << *this << endl
                << "\nFOAM aborting\n" << endl;
            ::abort();
        }
    }
}


Foam::Ostream& Foam::operator<<(Ostream& os, const error& fErr)
{
    os  << endl << fErr.message().c_str();

    if (error::level >= 2 && fErr.sourceFileLineNumber())
    {
        os  << endl << endl
            << "    From function " << fErr.functionName().c_str() << endl
            << "    in file " << fErr.sourceFileName().c_str()
            << " at line " << fErr.sourceFileLineNumber() << '.';
    }

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions

Foam::error Foam::FatalError("--> FOAM FATAL ERROR : ");

// ************************************************************************* //
