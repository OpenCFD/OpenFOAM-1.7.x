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
    Foam::basicThermoParcel

Description


SourceFiles
    basicThermoParcel.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermoParcel_H
#define basicThermoParcel_H

#include "ThermoParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class basicThermoParcel Declaration
\*---------------------------------------------------------------------------*/

class basicThermoParcel
:
    public ThermoParcel<basicThermoParcel>
{

public:

    //- Runtime type information
    TypeName("basicThermoParcel");

    // Constructors

       //- Construct from owner, position, and cloud owner
       //  Other properties initialised as null
       basicThermoParcel
       (
            ThermoCloud<basicThermoParcel>& owner,
            const vector position,
            const label cellI
       );

       //- Construct from components
       basicThermoParcel
       (
            ThermoCloud<basicThermoParcel>& owner,
            const vector position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector U0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        basicThermoParcel
        (
            const Cloud<basicThermoParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        basicThermoParcel(const basicThermoParcel& p);

        //- Construct and return a clone
        autoPtr<basicThermoParcel> clone() const
        {
            return
                autoPtr<basicThermoParcel>
                (
                    new basicThermoParcel(*this)
                );
        }


    //- Destructor
    virtual ~basicThermoParcel();
};


template<>
inline bool contiguous<basicThermoParcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
