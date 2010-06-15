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
    Foam::dispersionLESModel

Description
    Foam::dispersionLESModel

\*---------------------------------------------------------------------------*/

#ifndef dispersionLESModel_H
#define dispersionLESModel_H

#include "dispersionModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dispersionLESModel Declaration
\*---------------------------------------------------------------------------*/

class dispersionLESModel
:
    public dispersionModel
{
    // Private data

        const compressible::LESModel& turbulence_;


public:

    //- Runtime type information
    TypeName("dispersionLESModel");


    // Constructors

        //- Construct from components
        dispersionLESModel
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        virtual ~dispersionLESModel();


    // Member Functions

        const compressible::LESModel& turbulence() const
        {
            return turbulence_;
        }

        virtual void disperseParcels() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
