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

Class
    Foam::collisionModel

Description
    Base case for selecting a collision model

\*---------------------------------------------------------------------------*/

#ifndef collisionModel_H
#define collisionModel_H

#include "IOdictionary.H"
#include "spray.H"
#include "Random.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class collisionModel Declaration
\*---------------------------------------------------------------------------*/

class collisionModel
{

protected:

    // Protected data
        const dictionary& dict_;
        spray& spray_;
        Random& rndGen_;

        //dictionary coeffsDict_;


public:

    //- Runtime type information
    TypeName("collisionModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            collisionModel,
            dictionary,
            (
                const dictionary& dict,
                spray& sm,
                Random& rndGen
            ),
            (dict, sm, rndGen)
        );


    // Constructors

        //- Construct from components
        collisionModel
        (
            const dictionary& dict,
            spray& sm,
            Random& rndGen
        );


    // Destructor

        virtual ~collisionModel();


    // Selector

        static autoPtr<collisionModel> New
        (
            const dictionary& dict,
            spray& sm,
            Random& rndGen
        );


    // Member Functions

        virtual void collideParcels(const scalar dt) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
