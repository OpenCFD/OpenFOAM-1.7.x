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
    Foam::BasicReactingMultiphaseParcel

Description


SourceFiles
    BasicReactingMultiphaseParcel.C
    BasicReactingMultiphaseParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef BasicReactingMultiphaseParcel_H
#define BasicReactingMultiphaseParcel_H

#include "ReactingMultiphaseParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

// Forward declaration of classes
template<class ThermoType>
class BasicReactingMultiphaseParcel;

/*---------------------------------------------------------------------------*\
               Class BasicReactingMultiphaseParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class BasicReactingMultiphaseParcel
:
    public ReactingMultiphaseParcel<BasicReactingMultiphaseParcel<ThermoType> >
{

public:

    //- The type of thermodynamics this parcel was instantiated for
    typedef ThermoType thermoType;

    //- Run-time type information
    TypeName("BasicReactingMultiphaseParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        BasicReactingMultiphaseParcel
        (
             ReactingMultiphaseCloud<BasicReactingMultiphaseParcel>& owner,
             const vector& position,
             const label cellI
        );

        //- Construct from components
        BasicReactingMultiphaseParcel
        (
             ReactingMultiphaseCloud<BasicReactingMultiphaseParcel>& owner,
             const vector& position,
             const label cellI,
             const label typeId,
             const scalar nParticle0,
             const scalar d0,
             const vector& U0,
             const scalarField& YGas0,
             const scalarField& YLiquid0,
             const scalarField& YSolid0,
             const scalarField& Y0,
             const typename
                ReactingMultiphaseParcel<BasicReactingMultiphaseParcel>::
                constantProperties& constProps
        );

        //- Construct from Istream
        BasicReactingMultiphaseParcel
        (
            const Cloud<BasicReactingMultiphaseParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        BasicReactingMultiphaseParcel(const BasicReactingMultiphaseParcel& p);

        //- Construct and return a clone
        autoPtr<BasicReactingMultiphaseParcel> clone() const
        {
            return
                autoPtr<BasicReactingMultiphaseParcel>
                (
                    new BasicReactingMultiphaseParcel(*this)
                );
        }


    //- Destructor
    virtual ~BasicReactingMultiphaseParcel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicReactingMultiphaseParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
