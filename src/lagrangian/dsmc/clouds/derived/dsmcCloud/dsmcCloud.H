/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::dsmcCloud

Description
    Cloud class to simulate dsmc parcels

SourceFiles
    dsmcCloud.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCloud_H
#define dsmcCloud_H

#include "DsmcCloud.H"
#include "dsmcParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class dsmcCloud Declaration
\*---------------------------------------------------------------------------*/

class dsmcCloud
:
    public DsmcCloud<dsmcParcel>
{
    // Private member functions

        //- Disallow default bitwise copy construct
        dsmcCloud(const dsmcCloud&);

        //- Disallow default bitwise assignment
        void operator=(const dsmcCloud&);


public:

    //- Runtime type information
    TypeName("dsmcCloud");


    // Constructors

        //- Construct from components
        dsmcCloud
        (
            const word& cloudName,
            const volScalarField& T,
            const volVectorField& U
        );

        //- Construct from name and mesh, used to initialise.
        dsmcCloud
        (
            const word& cloudName,
            const fvMesh& mesh
        );

    //- Destructor
    ~dsmcCloud();


    // Member functions

        //- Write fields
        virtual void writeFields() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
