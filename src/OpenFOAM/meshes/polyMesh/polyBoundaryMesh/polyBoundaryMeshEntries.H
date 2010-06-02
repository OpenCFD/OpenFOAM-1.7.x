/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::polyBoundaryMeshEntries

Description
    Foam::polyBoundaryMeshEntries

SourceFiles
    polyBoundaryMeshEntries.C

\*---------------------------------------------------------------------------*/

#ifndef polyBoundaryMeshEntries_H
#define polyBoundaryMeshEntries_H

#include "regIOobject.H"
#include "PtrList.H"
#include "entry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyBoundaryMeshEntries
:
    public regIOobject,
    public PtrList<entry>
{

public:

    //- Runtime type information
    TypeName("polyBoundaryMesh");


    // Constructors

        explicit polyBoundaryMeshEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream(typeName))
        {
            close();
        }


   // Member functions

        bool writeData(Ostream&) const
        {
            notImplemented("writeData(Ostream&) const");
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
