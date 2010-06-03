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
    Foam::RanzMarshall

Description
    The Ranz-Marshal correlation for heat transfer

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshall_H
#define RanzMarshall_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class RanzMarshall Declaration
\*---------------------------------------------------------------------------*/

class RanzMarshall
:
    public heatTransferModel
{
private:

    // Private data

        dictionary heatDict_;
        scalar preRePrFactor_;
        scalar ReExponent_;
        scalar PrExponent_;

public:

    //- Runtime type information
        TypeName("RanzMarshall");


    // Constructors

        //- Construct from dictionary
        RanzMarshall
        (
            const dictionary& dict
        );


    // Destructor

        ~RanzMarshall();


    // Member Functions

        bool heatTransfer() const;

        //- Nusselt Number
        scalar Nu
        (
            const scalar ReynoldsNumber,
            const scalar PrandtlNumber
        ) const;

        scalar relaxationTime
        (
            const scalar liquidDensity,
            const scalar diameter,
            const scalar liquidcL,
            const scalar kappa,
            const scalar ReynoldsNumber,
            const scalar PrandtlNumber
        ) const;

        scalar fCorrection(const scalar z) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
