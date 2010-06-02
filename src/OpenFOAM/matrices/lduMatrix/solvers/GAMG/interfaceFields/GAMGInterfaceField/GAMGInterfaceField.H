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
    Foam::GAMGInterfaceField

Description
     Abstract base class for GAMG agglomerated interface fields.

SourceFiles
    GAMGInterfaceField.C
    newAmgInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGInterfaceField_H
#define GAMGInterfaceField_H

#include "lduInterfaceField.H"
#include "GAMGInterface.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class GAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class GAMGInterfaceField
:
    public lduInterfaceField
{
    // Private data

        //- Local reference cast into the interface
        const GAMGInterface& interface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        GAMGInterfaceField(const GAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const GAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("GAMGInterfaceField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGInterfaceField,
            lduInterface,
            (
                const GAMGInterface& GAMGCp,
                const lduInterfaceField& fineInterface
            ),
            (GAMGCp, fineInterface)
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<GAMGInterfaceField> New
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterface
        );


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        GAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField&
        )
        :
            lduInterfaceField(GAMGCp),
            interface_(GAMGCp)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
