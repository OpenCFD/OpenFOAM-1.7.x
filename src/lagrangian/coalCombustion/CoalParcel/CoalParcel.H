/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    CoalParcel

Description


SourceFiles
    CoalParcel.C
    CoalParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef CoalParcel_H
#define CoalParcel_H

#include "ReactingMultiphaseParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class ThermoType>
class CoalParcel;

/*---------------------------------------------------------------------------*\
                         Class CoalParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class CoalParcel
:
    public ReactingMultiphaseParcel<CoalParcel<ThermoType> >
{
public:

    //- The type of thermodynamics this parcel was instantiated for
    typedef ThermoType thermoType;

    // Run-time type information
    TypeName("CoalParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        CoalParcel
        (
             ReactingMultiphaseCloud<CoalParcel>& owner,
             const vector& position,
             const label cellI
        );

        //- Construct from components
        CoalParcel
        (
             ReactingMultiphaseCloud<CoalParcel>& owner,
             const vector& position,
             const label cellI,
             const label typeId,
             const scalar nParticle0,
             const scalar d0,
             const vector& U0,
             const scalarField& YMixture0,
             const scalarField& YGas0,
             const scalarField& YLiquid0,
             const scalarField& YSolid0,
             const typename
                 ReactingMultiphaseParcel<CoalParcel>::
                 constantProperties& constProps
        );

        //- Construct from Istream
        CoalParcel
        (
            const Cloud<CoalParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<CoalParcel> clone() const
        {
            return autoPtr<CoalParcel>(new CoalParcel(*this));
        }


    //- Destructor
    virtual ~CoalParcel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CoalParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
