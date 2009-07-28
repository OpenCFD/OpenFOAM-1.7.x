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
    Foam::RanzMarshall

Description
    The Ranz-Marshall correlation for heat transfer

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshall_H
#define RanzMarshall_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class RanzMarshall Declaration
\*---------------------------------------------------------------------------*/

template <class CloudType>
class RanzMarshall
:
    public HeatTransferModel<CloudType>
{
    // Private data

        // Prandtl number
        const scalar Pr_;


public:

    //- Runtime type information
    TypeName("RanzMarshall");


    // Constructors

        //- Construct from dictionary
        RanzMarshall
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~RanzMarshall();


    // Member Functions

        //- Flag to indicate whether model activates heat transfer model
        virtual bool active() const;

        //- Nusselt number
        virtual scalar Nu
        (
            const scalar Re,
            const scalar Pr
        ) const;

        //- Prandtl number
        virtual scalar Pr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "RanzMarshall.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
