/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::LiquidEvaporation

Description
    Liquid evaporation model
    - uses ideal gas assumption

\*---------------------------------------------------------------------------*/

#ifndef LiquidEvaporation_H
#define LiquidEvaporation_H

#include "PhaseChangeModel.H"
#include "liquidMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class LiquidEvaporation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LiquidEvaporation
:
    public PhaseChangeModel<CloudType>
{
protected:

    // Protected data

        //- Global liquid properties data
        autoPtr<liquidMixture> liquids_;

        //- List of active liquid names
        List<word> activeLiquids_;

        //- Mapping between liquid and carrier species
        List<label> liqToCarrierMap_;

        //- Mapping between local and global liquid species
        List<label> liqToLiqMap_;


    // Protected member functions

        //- Sherwood number as a function of Reynolds and Schmidt numbers
        scalar Sh(const scalar Re, const scalar Sc) const;

        //- Calculate the carrier phase component volume fractions at cellI
        scalarField calcXc(const label cellI) const;


public:

    //- Runtime type information
    TypeName("LiquidEvaporation");


    // Constructors

        //- Construct from dictionary
        LiquidEvaporation
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~LiquidEvaporation();


    // Member Functions

        //- Flag to indicate whether model activates phase change model
        virtual bool active() const;

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const label cellI,
            const scalar Re,
            const scalar d,
            const scalar nu,
            const scalar T,
            const scalar Ts,
            const scalar pc,
            scalarField& dMassPC
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LiquidEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
