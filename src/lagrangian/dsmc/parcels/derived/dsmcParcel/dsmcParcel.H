/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::dsmcParcel

Description


SourceFiles
    dsmcParcel.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcParcel_H
#define dsmcParcel_H

#include "DsmcParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dsmcParcel Declaration
\*---------------------------------------------------------------------------*/

class dsmcParcel
:
    public DsmcParcel<dsmcParcel>
{

public:

    //- Run-time type information
    TypeName("dsmcParcel");


    // Constructors

        //- Construct from components
        dsmcParcel
        (
            DsmcCloud<dsmcParcel>& owner,
            const vector& position,
            const vector& U,
            const scalar Ei,
            const label celli,
            const label typeId
        );

        //- Construct from Istream
        dsmcParcel
        (
            const Cloud<dsmcParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<dsmcParcel> clone() const
        {
            return autoPtr<dsmcParcel>(new dsmcParcel(*this));
        }


    //- Destructor
    virtual ~dsmcParcel();
};


template<>
inline bool contiguous<dsmcParcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
