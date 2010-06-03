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

#include "regIOobject.H"
#include "IFstream.H"
#include "Time.H"
#include "PstreamReduceOps.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::regIOobject::readStream()
{
    if (IFstream::debug)
    {
        Info<< "regIOobject::readStream() : "
            << "reading object " << name()
            << " from file " << objectPath()
            << endl;
    }

    if (readOpt() == NO_READ)
    {
        FatalErrorIn("regIOobject::readStream()")
            << "NO_READ specified for read-constructor of object " << name()
            << " of class " << headerClassName()
            << abort(FatalError);
    }

    // Construct object stream and read header if not already constructed
    if (!isPtr_)
    {
        if (!(isPtr_ = objectStream()))
        {
            FatalIOError
            (
                "regIOobject::readStream()",
                __FILE__,
                __LINE__,
                objectPath(),
                0
            )   << "cannot open file"
                << exit(FatalIOError);
        }
        else if (!readHeader(*isPtr_))
        {
            FatalIOErrorIn("regIOobject::readStream()", *isPtr_)
                << "problem while reading header for object " << name()
                << exit(FatalIOError);
        }
    }

    if (!lastModified_)
    {
        lastModified_ = lastModified(filePath());
    }

    return *isPtr_;
}


Foam::Istream& Foam::regIOobject::readStream(const word& expectName)
{
    if (IFstream::debug)
    {
        Info<< "regIOobject::readStream(const word&) : "
            << "reading object " << name()
            << " from file " << objectPath()
            << endl;
    }

    // Construct IFstream if not already constructed
    if (!isPtr_)
    {
        readStream();

        // Check the className of the regIOobject
        // dictionary is an allowable name in case the actual class
        // instantiated is a dictionary
        if
        (
            expectName.size()
         && headerClassName() != expectName
         && headerClassName() != "dictionary"
        )
        {
            FatalIOErrorIn("regIOobject::readStream(const word&)", *isPtr_)
                << "unexpected class name " << headerClassName()
                << " expected " << expectName << endl
                << "    while reading object " << name()
                << exit(FatalIOError);
        }
    }

    return *isPtr_;
}


void Foam::regIOobject::close()
{
    if (IFstream::debug)
    {
        Info<< "regIOobject::close() : "
            << "finished reading " << filePath()
            << endl;
    }

    if (isPtr_)
    {
        delete isPtr_;
        isPtr_ = NULL;
    }
}


bool Foam::regIOobject::readData(Istream&)
{
    return false;
}


bool Foam::regIOobject::read()
{
    bool ok = readData(readStream(type()));
    close();
    return ok;
}


bool Foam::regIOobject::modified() const
{
    return
    (
        lastModified_
     && lastModified(filePath()) > (lastModified_ + fileModificationSkew)
    );
}


bool Foam::regIOobject::readIfModified()
{
    if (lastModified_)
    {
        time_t newTimeStamp = lastModified(filePath());

        bool readFile = false;

        if (newTimeStamp > (lastModified_ + fileModificationSkew))
        {
            readFile = true;
        }

        if (Pstream::parRun())
        {
            bool readFileOnThisProc = readFile;
            reduce(readFile, andOp<bool>());

            if (readFileOnThisProc && !readFile)
            {
                WarningIn("regIOobject::readIfModified()")
                    << "Delaying reading " << name()
                    << " of class " << headerClassName()
                    << " due to inconsistent "
                       "file time-stamps between processors"
                    << endl;
            }
        }

        if (readFile)
        {
            lastModified_ = newTimeStamp;
            Info<< "regIOobject::readIfModified() : " << nl
                << "    Reading object " << name()
                << " from file " << filePath() << endl;
            return read();
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
