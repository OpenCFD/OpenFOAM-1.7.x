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

#include "dynamicFvMesh.H"
#include "Time.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicFvMesh> Foam::dynamicFvMesh::New(const IOobject& io)
{
    // Enclose the creation of the dynamicMesh to ensure it is
    // deleted before the dynamicFvMesh is created otherwise the dictionary
    // is entered in the database twice
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            io.time().constant(),
            io.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    word dynamicFvMeshTypeName(dynamicMeshDict.lookup("dynamicFvMesh"));

    Info<< "Selecting dynamicFvMesh " << dynamicFvMeshTypeName << endl;

    dlLibraryTable::open
    (
        dynamicMeshDict,
        "dynamicFvMeshLibs",
        IOobjectConstructorTablePtr_
    );

    if (!IOobjectConstructorTablePtr_)
    {
        FatalErrorIn
        (
            "dynamicFvMesh::New(const IOobject&)"
        )   << "dynamicFvMesh table is empty"
            << exit(FatalError);
    }

    IOobjectConstructorTable::iterator cstrIter =
        IOobjectConstructorTablePtr_->find(dynamicFvMeshTypeName);

    if (cstrIter == IOobjectConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dynamicFvMesh::New(const IOobject&)"
        )   << "Unknown dynamicFvMesh type " << dynamicFvMeshTypeName
            << endl << endl
            << "Valid dynamicFvMesh types are :" << endl
            << IOobjectConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<dynamicFvMesh>(cstrIter()(io));
}


// ************************************************************************* //
