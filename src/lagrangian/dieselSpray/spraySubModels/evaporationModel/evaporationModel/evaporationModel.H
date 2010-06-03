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
    Foam::evaporationModel

Description
    Base class for selecting an evaporation model

\*---------------------------------------------------------------------------*/

#ifndef evaporationModel_H
#define evaporationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class evaporationModel Declaration
\*---------------------------------------------------------------------------*/

class evaporationModel
{
protected:

    // Protected data

        const dictionary& dict_;

public:

    //- Runtime type information
        TypeName("evaporationModel");

    //- Evaluate the gaseous properties Tand Y (implicit/explicit)
    //word gasPropertyEvaluation_;
    //label nEvapIter_;


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            evaporationModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from components
        evaporationModel
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~evaporationModel();


    // Selector

        static autoPtr<evaporationModel> New
        (
            const dictionary& dict
        );


    // Member Functions

    //- evaporation on/off
        virtual bool evaporation() const = 0;

        //- The Sherwood number
        virtual scalar Sh
        (
            const scalar ReynoldsNumber,
            const scalar SchmidtNumber
        ) const = 0;

        virtual scalar relaxationTime
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
        ) const = 0;

        virtual scalar boilingTime
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
            const scalar kLiq
        ) const = 0;

        virtual label nEvapIter() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
