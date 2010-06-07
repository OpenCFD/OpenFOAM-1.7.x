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
    Foam::standardEvaporationModel

Description
    standard evaporation model based on the logarithmic expression
    of the ration of pressure differences.

\*---------------------------------------------------------------------------*/

#ifndef standardEvaporationModel_H
#define standardEvaporationModel_H

#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class standardEvaporationModel Declaration
\*---------------------------------------------------------------------------*/

class standardEvaporationModel
:
    public evaporationModel
{

private:

    // Private data

        dictionary evapDict_;
        scalar preReScFactor_;
        scalar ReExponent_;
        scalar ScExponent_;
    
        word evaporationScheme_;
        label nEvapIter_;


public:

    //- Runtime type information
        TypeName("standardEvaporationModel");


    // Constructors

        //- Construct from dictionary
        standardEvaporationModel
        (
            const dictionary& dict
        );


    // Destructor

        ~standardEvaporationModel();


    // Member Functions

        bool evaporation() const;

       //- Correlation for the Sherwood Number
        scalar Sh
        (
            const scalar ReynoldsNumber,
            const scalar SchmidtNumber
        ) const;

        //- Return the evaporation relaxation time
        scalar relaxationTime
        (
            const scalar diameter,
            const scalar liquidDensity,
            const scalar rhoFuelVapor,
            const scalar massDiffusionCoefficient,
            const scalar ReynoldsNumber,
            const scalar SchmidtNumber,
            const scalar Xs,
            const scalar Xf,
            const scalar m0,
            const scalar dm,
            const scalar dt
        ) const;

        scalar boilingTime
        (
            const scalar liquidDensity,
            const scalar cpFuel,
            const scalar heatOfVapour,
            const scalar kappa,
            const scalar Nusselt,
            const scalar deltaTemp,
            const scalar diameter,
            const scalar, 
            const scalar, 
            const scalar, 
            const scalar, 
            const scalar, 
            const scalar, 
            const scalar, 
            const scalar, 
            const scalar 
        ) const;

        inline label nEvapIter() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
