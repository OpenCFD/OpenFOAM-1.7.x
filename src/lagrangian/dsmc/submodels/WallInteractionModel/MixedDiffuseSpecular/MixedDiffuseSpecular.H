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
    Foam::MixedDiffuseSpecular

Description
    Wall interaction setting microscopic velocity to a random one drawn from a
    Maxwellian distribution corresponding to a specified temperature

\*---------------------------------------------------------------------------*/

#ifndef MixedDiffuseSpecular_H
#define MixedDiffuseSpecular_H

#include "WallInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class MixedDiffuseSpecular Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MixedDiffuseSpecular
:
    public WallInteractionModel<CloudType>
{
    // Private data

        //- Fraction of wall interactions that are diffuse
        scalar diffuseFraction_;


public:

    //- Runtime type information
    TypeName("MixedDiffuseSpecular");


    // Constructors

        //- Construct from dictionary
        MixedDiffuseSpecular
        (
            const dictionary& dict,
            CloudType& cloud
        );


    // Destructor
    virtual ~MixedDiffuseSpecular();


    // Member Functions

        //- Apply wall correction
        virtual void correct
        (
            const wallPolyPatch& wpp,
            const label faceId,
            vector& U,
            scalar& Ei,
            label typeId
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MixedDiffuseSpecular.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
