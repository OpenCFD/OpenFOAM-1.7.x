/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::basicKinematicParcel

Description


SourceFiles
    basicKinematicParcel.C

\*---------------------------------------------------------------------------*/

#ifndef basicKinematicParcel_H
#define basicKinematicParcel_H

#include "KinematicParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class basicKinematicParcel Declaration
\*---------------------------------------------------------------------------*/

class basicKinematicParcel
:
    public KinematicParcel<basicKinematicParcel>
{

public:

    //- Run-time type information
    TypeName("basicKinematicParcel");


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        basicKinematicParcel
        (
            KinematicCloud<basicKinematicParcel>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        basicKinematicParcel
        (
            KinematicCloud<basicKinematicParcel>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        basicKinematicParcel
        (
            const Cloud<basicKinematicParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        basicKinematicParcel(const basicKinematicParcel& p);

        //- Construct and return a clone
        autoPtr<basicKinematicParcel> clone() const
        {
            return
                autoPtr<basicKinematicParcel>
                (
                    new basicKinematicParcel(*this)
                );
        }


    //- Destructor
    virtual ~basicKinematicParcel();
};


template<>
inline bool contiguous<basicKinematicParcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
