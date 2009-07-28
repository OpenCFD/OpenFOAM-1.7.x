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
    Foam::AverageIOField

Description
    A primitive field + average with IO.

SourceFiles
    AverageIOField.C

\*---------------------------------------------------------------------------*/

#ifndef AverageIOField_H
#define AverageIOField_H

#include "regIOobject.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AverageIOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class AverageIOField
:
    public regIOobject,
    public pTraits<Type>,
    public Field<Type>
{

public:

    TypeName("AverageField");


    // Constructors

        //- Construct from IOobject
        AverageIOField
        (
            const IOobject&
        );

        //- Construct from IOobject and size (does not set values)
        AverageIOField
        (
            const IOobject&,
            const label size
        );

        //- Construct from components
        AverageIOField
        (
            const IOobject&,
            const Type& average,
            const Field<Type>&
        );


    // Member functions

        const pTraits<Type>& average() const
        {
            return static_cast<const pTraits<Type>&>(*this);
        }

        bool writeData(Ostream& os) const
        {
            os  << static_cast<const Type&>(*this)
                << token::NL
                << static_cast<const Field<Type>&>(*this);

            return os.good();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "AverageIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
