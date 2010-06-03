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
    Foam::RutlandFlashBoil

Description
    Evaporation model with flash boiling.

    Evaporation model that takes into account of:
      - flash boiling of the liquid core and droplets
      - superheated vaporization of a boiling droplet

    Accurated description in
    @verbatim
    B. Zuo, A.M. Gomes, C.J. Rutland
    "Modeling Superheated Fuel Spray Vaporization"
    Int. Journal of Engine Research, 2000. Vol. 1, pp. 321-326
    @endverbatim
\*---------------------------------------------------------------------------*/

#ifndef RutlandFlashBoil_H
#define RutlandFlashBoil_H

#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class RutlandFlashBoil Declaration
\*---------------------------------------------------------------------------*/

class RutlandFlashBoil
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
        TypeName("RutlandFlashBoil");


    // Constructors

        //- Construct from dictionary
        RutlandFlashBoil
        (
            const dictionary& dict
        );


    // Destructor

        ~RutlandFlashBoil();


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
            const scalar liquidCore,
            const scalar time,
            const scalar tDrop,
            const scalar tBoilingSurface,
            const scalar vapourSurfaceEnthalpy,
            const scalar vapourFarEnthalpy,
            const scalar cpGas,
            const scalar temperature,
            const scalar kLiquid
        ) const;

        inline label nEvapIter() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
