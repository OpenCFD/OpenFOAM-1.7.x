/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::BasicReactingParcel

Description


SourceFiles
    BasicReactingParcel.C
    BasicReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef BasicReactingParcel_H
#define BasicReactingParcel_H

#include "ReactingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class ThermoType>
class BasicReactingParcel;

/*---------------------------------------------------------------------------*\
                    Class BasicReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class BasicReactingParcel
:
    public ReactingParcel<BasicReactingParcel<ThermoType> >
{

public:

    //- The type of thermodynamics this parcel was instantiated for
    typedef ThermoType thermoType;

    //- Run-time type information
    TypeName("BasicReactingParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        BasicReactingParcel
        (
            ReactingCloud<BasicReactingParcel>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        BasicReactingParcel
        (
            ReactingCloud<BasicReactingParcel>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const scalarField& Y0,
            const typename ReactingParcel<BasicReactingParcel>::
                constantProperties& constProps
        );

        //- Construct from Istream
        BasicReactingParcel
        (
            const Cloud<BasicReactingParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        BasicReactingParcel(const BasicReactingParcel& p);

        //- Construct and return a clone
        autoPtr<BasicReactingParcel> clone() const
        {
            return
                autoPtr<BasicReactingParcel>
                (
                    new BasicReactingParcel(*this)
                );
        }


    //- Destructor
    virtual ~BasicReactingParcel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
