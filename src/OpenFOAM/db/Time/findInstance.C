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
    If "name" is empty: return the location of "directory"
    If "name" is not empty: return the location of "directory" containing the
    file "name".
    Used in reading mesh data.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::Time::findInstance
(
    const fileName& dir,
    const word& name,
    const IOobject::readOption rOpt
) const
{
    // Note: if name is empty, just check the directory itself

    // check the current time directory
    if
    (
        name.empty()
      ? isDir(path()/timeName()/dir)
      :
        (
            isFile(path()/timeName()/dir/name)
         && IOobject(name, timeName(), dir, *this).headerOk()
        )
    )
    {
        if (debug)
        {
            Info<< "Time::findInstance"
                "(const fileName&, const word&, const IOobject::readOption)"
                << " : found \"" << name
                << "\" in " << timeName()/dir
                << endl;
        }

        return timeName();
    }

    // Search back through the time directories to find the time
    // closest to and lower than current time

    instantList ts = times();
    label instanceI;

    for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
    {
        if (ts[instanceI].value() <= timeOutputValue())
        {
            break;
        }
    }

    // continue searching from here
    for (; instanceI >= 0; --instanceI)
    {
        if
        (
            name.empty()
          ? isDir(path()/ts[instanceI].name()/dir)
          :
            (
                isFile(path()/ts[instanceI].name()/dir/name)
             && IOobject(name, ts[instanceI].name(), dir, *this).headerOk()
            )
        )
        {
            if (debug)
            {
                Info<< "Time::findInstance"
                    "(const fileName&, const word&, const IOobject::readOption)"
                    << " : found \"" << name
                    << "\" in " << ts[instanceI].name()/dir
                    << endl;
            }

            return ts[instanceI].name();
        }
    }


    // not in any of the time directories, try constant

    // Note. This needs to be a hard-coded constant, rather than the
    // constant function of the time, because the latter points to
    // the case constant directory in parallel cases

    if
    (
        name.empty()
      ? isDir(path()/constant()/dir)
      :
        (
            isFile(path()/constant()/dir/name)
         && IOobject(name, constant(), dir, *this).headerOk()
        )
    )
    {
        if (debug)
        {
            Info<< "Time::findInstance"
                "(const fileName&, const word&, const IOobject::readOption)"
                << " : found \"" << name
                << "\" in " << constant()/dir
                << endl;
        }

        return constant();
    }

    if (rOpt == IOobject::MUST_READ)
    {
        FatalErrorIn
        (
            "Time::findInstance"
            "(const fileName&, const word&, const IOobject::readOption)"
        )   << "Cannot find file \"" << name << "\" in directory "
            << constant()/dir
            << exit(FatalError);
    }

    return constant();
}


// ************************************************************************* //
