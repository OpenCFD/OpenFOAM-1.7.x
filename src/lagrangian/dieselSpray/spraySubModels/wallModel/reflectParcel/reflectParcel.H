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
    Foam::reflectParcel

Description
    The parcels are reflected back from the wall in an inelastic collision.

\*---------------------------------------------------------------------------*/

#ifndef reflectParcel_H
#define reflectParcel_H

#include "wallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reflectParcel Declaration
\*---------------------------------------------------------------------------*/

class reflectParcel
:
    public wallModel
{

private:
    // Private data

        const volVectorField& U_;

        dictionary coeffsDict_;

    // model constants

    //- Elasticity constant v^2_out = elasticity * v^2_in
        scalar elasticity_;

public:

    //- Runtime type information
    TypeName("reflect");


    // Constructors

        //- Construct from components
        reflectParcel
        (
            const dictionary& dict,
            const volVectorField& U,
            spray& sm
        );


    // Destructor

        virtual ~reflectParcel();


    // Member Functions

        //- Return true if parcel is to be kept, and false if it is to be
        // removed
        bool wallTreatment
        (
            parcel& parcel,
            const label facei
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
