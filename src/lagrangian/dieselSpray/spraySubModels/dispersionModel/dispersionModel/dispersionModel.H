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
    Foam::dispersionModel

Description
    Foam::dispersionModel

\*---------------------------------------------------------------------------*/

#ifndef dispersionModel_H
#define dispersionModel_H

#include "IOdictionary.H"
#include "spray.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dispersionModel Declaration
\*---------------------------------------------------------------------------*/

class dispersionModel
{

protected:

    // Protected data

        const dictionary& dict_;
        spray& spray_;


public:

    //- Runtime type information
    TypeName("dispersionModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dispersionModel,
            dictionary,
            (
                const dictionary& dict,
                spray& sm
            ),
            (dict, sm)
        );


    // Constructors

        //- Construct from components
        dispersionModel
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        virtual ~dispersionModel();


    // Selector

        static autoPtr<dispersionModel> New
        (
            const dictionary& dict,
            spray& sm
        );


    // Member Functions

        virtual void disperseParcels() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
