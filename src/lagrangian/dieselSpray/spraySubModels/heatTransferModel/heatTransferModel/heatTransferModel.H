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
    Foam::heatTransferModel

Description
    Foam::heatTransferModel

\*---------------------------------------------------------------------------*/

#ifndef heatTransferModel_H
#define heatTransferModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferModel

{

protected:

    // Protected data

    const dictionary& dict_;

public:

    //- Runtime type information
        TypeName("heatTransferModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from dictionary
        heatTransferModel
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~heatTransferModel();


    // Selector

        static autoPtr<heatTransferModel> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- heat transfer on/off
        virtual bool heatTransfer() const = 0;

        //- Return Nusselt Number
        virtual scalar Nu
        (
            const scalar ReynoldsNumber,
            const scalar PrandtlNumber
        ) const = 0;

        //- Return the heat transfer relaxation time
        virtual scalar relaxationTime
        (
            const scalar liquidDensity,
            const scalar diameter,
            const scalar liquidcL,
            const scalar kappa,
            const scalar ReynoldsNumber,
            const scalar PrandtlNumber
        ) const = 0;

        //- Return the correction function f...
        virtual scalar fCorrection(const scalar z) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
