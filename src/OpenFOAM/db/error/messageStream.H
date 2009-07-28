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
    Foam::messageStream

Description
    Class to handle messaging in a simple, consistent stream-based
    manner.

    The messageStream class is globaly instantiated with a title string a
    given severity, which controls the program termination, and a number of
    errors before termination.  Errors, messages and other data are piped to
    the messageStream class in the standard manner.

Usage
    @code
        messageStream << "message1" << "message2" << FoamDataType << endl;
    @endcode

SourceFiles
    messageStream.C

\*---------------------------------------------------------------------------*/

#ifndef messageStream_H
#define messageStream_H

#include "label.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOstream;
class Ostream;
class OSstream;
class OStringStream;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class messageStream Declaration
\*---------------------------------------------------------------------------*/

class messageStream
{

public:

    //- Severity flags
    enum errorSeverity
    {
        INFO,       // Debugging information in event of error
        WARNING,    // Warning of possible problem
        SERIOUS,    // A serious problem (data corruption?)
        FATAL       // Oh bugger!
    };


protected:

    // Private data

        string title_;
        errorSeverity severity_;
        int maxErrors_;
        int errorCount_;


public:

    // Debug switches

        static int level;


    // Constructors

        //- Construct from components
        messageStream
        (
            const string& title,
            errorSeverity sev,
            const int maxErrors = 0
        );


        //- Construct from dictionary
        messageStream(const dictionary& dict);


    // Member functions

        //- Return the title of this error type
        const string& title() const
        {
            return title_;
        }

        //- Return the maximum number of errors before program termination
        int maxErrors() const
        {
            return maxErrors_;
        }

        //- Return non-const access to the maximum number of errors before
        //  program termination to enable user to reset it
        int& maxErrors()
        {
            return maxErrors_;
        }

        //- Convert to Ostream
        //  Prints basic message and then returns Ostream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        OSstream& operator()
        (
            const string& functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Convert to Ostream
        //  Prints basic message and then returns Ostream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const string& ioFileName,
            const label ioStartLineNumber = -1,
            const label ioEndLineNumber = -1
        );

        //- Convert to Ostream
        //  Prints basic message and then returns Ostream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const IOstream&
        );

        //- Convert to Ostream
        //  Prints basic message and then returns Ostream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const dictionary&
        );

        //- Convert to Ostream for << operations
        operator OSstream&();

        //- Explicitly convert to Ostream for << operations
        OSstream& operator()()
        {
            return operator OSstream&();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error declarations: defined in messageStream.C

extern messageStream SeriousError;
extern messageStream Warning;
extern messageStream Info;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Convienient macros to add the file name and line number to the function name

#define SeriousErrorIn(fn) SeriousError(fn, __FILE__, __LINE__)
#define SeriousIOErrorIn(fn, ios) SeriousError(fn, __FILE__, __LINE__, ios)

#define WarningIn(fn) Warning(fn, __FILE__, __LINE__)
#define IOWarningIn(fn, ios) Warning(fn, __FILE__, __LINE__, ios)

#define InfoIn(fn) Info(fn, __FILE__, __LINE__)
#define IOInfoIn(fn, ios) Info(fn, __FILE__, __LINE__, ios)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
