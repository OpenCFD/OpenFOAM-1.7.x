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
    Foam::injectorModel


Description
    The injectormodel is used to set initial droplet size and the
    velocity direction 

\*---------------------------------------------------------------------------*/

#ifndef injectorModel_H
#define injectorModel_H

#include "liquid.H"
#include "injector.H"
#include "PtrList.H"
#include "IOdictionary.H"
#include "volFieldsFwd.H"
#include "autoPtr.H"
#include "spray.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class injectorModel Declaration
\*---------------------------------------------------------------------------*/

class injectorModel
{

protected:

    // Protected data

        const dictionary& dict_;
        spray& sm_;

        const PtrList<injector>& injectors_;    
        Random& rndGen_;

public:

    //-Runtime type information
        TypeName("injectorModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            injectorModel,
            dictionary,
            (
                const dictionary& dict,
                spray& sm
            ),
            (dict, sm)
        );


    // Constructors

        //- Construct from components
        injectorModel
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        virtual ~injectorModel();


    // Selector

        static autoPtr<injectorModel> New
        (
            const dictionary& dict,
            spray& sm
        );


    // Member Functions

        //- Return the initial droplet diameter
        virtual scalar d0(const label injector, const scalar time) const = 0;

        //- Return the spray angle of the injector
        virtual vector direction
        (
            const label injector, 
            const label hole,
            const scalar time,
            const scalar d
        ) const = 0;

        virtual scalar velocity
        (
            const label i,
            const scalar time
        ) const = 0;
        
        virtual scalar averageVelocity
        (
            const label i
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
