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
    regionProperties

Description
    Simple class to hold region information for coupled region simulations

SourceFiles
    regionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef regionProperties_H
#define regionProperties_H

#include "IOdictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionProperties Declaration
\*---------------------------------------------------------------------------*/

class regionProperties
:
    public IOdictionary
{
    // Private data

        //- List of the fluid region names
        List<word> fluidRegionNames_;

        //- List of the solid region names
        List<word> solidRegionNames_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionProperties(const regionProperties&);

        //- Disallow default bitwise assignment
        void operator=(const regionProperties&);


public:

    // Constructors

        //- Construct from components
        regionProperties(const Time& runTime);


    // Destructor

        ~regionProperties();


    // Member Functions

        // Access

            //- Return const reference to the list of fluid region names
            const List<word>& fluidRegionNames() const;

            //- Return const reference to the list of solid region names
            const List<word>& solidRegionNames() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
