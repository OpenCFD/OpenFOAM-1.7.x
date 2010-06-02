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
    Foam::ash

Description
    Coal ash

SourceFiles
    ash.C

\*---------------------------------------------------------------------------*/

#ifndef solid_ash_H
#define solid_ash_H

#include "solid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class ash;

Ostream& operator<<
(
     Ostream&,
     const ash&
);


/*---------------------------------------------------------------------------*\
                              Class ash Declaration
\*---------------------------------------------------------------------------*/

class ash
:
    public solid
{

public:

    //- Runtime type information
    TypeName("ash");


    // Constructors

        //- Construct null
        ash();

        //- Construct from solid
        ash(const solid& s);

        //- Construct from Istream
        ash(Istream& is);


    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const ash& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
