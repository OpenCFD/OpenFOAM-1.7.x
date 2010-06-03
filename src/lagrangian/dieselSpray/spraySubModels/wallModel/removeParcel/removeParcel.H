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
    Foam::removeParcel

Description
    remove parcel, regardless of boundary type.

\*---------------------------------------------------------------------------*/

#ifndef removeParcel_H
#define removeParcel_H

#include "wallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class removeParcel Declaration
\*---------------------------------------------------------------------------*/

class removeParcel
:
    public wallModel
{

public:

    //- Runtime type information
    TypeName("remove");


    // Constructors

        //- Construct from components
        removeParcel
        (
            const dictionary& dict,
            const volVectorField& U,
            spray& sm
        );


    // Destructor

        virtual ~removeParcel();


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
