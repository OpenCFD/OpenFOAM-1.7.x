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

#include "IOobjectList.H"
#include "Time.H"
#include "OSspecific.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobjectList::IOobjectList(const label nIoObjects)
:
    HashPtrTable<IOobject>(nIoObjects)
{}


Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local
)
:
    HashPtrTable<IOobject>()
{
    word newInstance = instance;

    if (!isDir(db.path(instance)))
    {
        newInstance = db.time().findInstancePath(instant(instance));

        if (newInstance.empty())
        {
            return;
        }
    }

    // Create list file names in directory
    fileNameList ObjectNames =
        readDir(db.path(newInstance, db.dbDir()/local), fileName::FILE);

    forAll(ObjectNames, i)
    {
        IOobject* objectPtr = new IOobject
        (
            ObjectNames[i],
            newInstance,
            local,
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (objectPtr->headerOk())
        {
            insert(ObjectNames[i], objectPtr);
        }
        else
        {
            delete objectPtr;
        }
    }
}


Foam::IOobjectList::IOobjectList(const IOobjectList& ioOL)
:
    HashPtrTable<IOobject>(ioOL)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::IOobjectList::~IOobjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOobjectList::add(IOobject& io)
{
    return insert(io.name(), &io);
}


bool Foam::IOobjectList::remove(IOobject& io)
{
    HashPtrTable<IOobject>::iterator iter =
        HashPtrTable<IOobject>::find(io.name());

    if (iter != end())
    {
        return erase(iter);
    }
    else
    {
        return false;
    }
}


Foam::IOobject* Foam::IOobjectList::lookup(const word& name) const
{
    HashPtrTable<IOobject>::const_iterator iter = find(name);

    if (iter != end())
    {
        if (IOobject::debug)
        {
            Info<< "IOobjectList::lookup : found " << name
                << endl;
        }

        return const_cast<IOobject*>(*iter);
    }
    else
    {
        if (IOobject::debug)
        {
            Info<< "IOobjectList::lookup : could not find " << name
                << endl;
        }

        return NULL;
    }
}


Foam::IOobjectList Foam::IOobjectList::lookupClass(const word& ClassName) const
{
    IOobjectList IOobjectsOfClass(size());

    for
    (
        HashPtrTable<IOobject>::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        if (iter()->headerClassName() == ClassName)
        {
            if (IOobject::debug)
            {
                Info<< "IOobjectList::lookupClass : found "
                    << iter()->name()
                    << endl;
            }

            IOobjectsOfClass.insert(iter()->name(), new IOobject(*iter()));
        }
    }

    return IOobjectsOfClass;
}


Foam::wordList Foam::IOobjectList::names() const
{
    wordList objectNames(size());

    label count = 0;
    for
    (
        HashPtrTable<IOobject>::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        objectNames[count++] = iter()->name();
    }

    return objectNames;
}


Foam::wordList Foam::IOobjectList::names(const word& ClassName) const
{
    wordList objectNames(size());

    label count = 0;
    for
    (
        HashPtrTable<IOobject>::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        if (iter()->headerClassName() == ClassName)
        {
            objectNames[count++] = iter()->name();
        }
    }

    objectNames.setSize(count);

    return objectNames;
}


// ************************************************************************* //
