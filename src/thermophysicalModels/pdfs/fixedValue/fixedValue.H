/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::fixedValue

Description
    Returns a fixed value

SourceFiles
    fixedValue.C

\*---------------------------------------------------------------------------*/

#ifndef pdfFixedValue_H
#define pdfFixedValue_H

#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pdfs
{
/*---------------------------------------------------------------------------*\
                       Class fixedValue Declaration
\*---------------------------------------------------------------------------*/

class fixedValue
:
    public pdf
{
    // Private data

        //- Fixed value
        scalar value_;


public:

    //- Runtime type information
    TypeName("fixedValue");


    // Constructors

        //- Construct from components
        fixedValue
        (
            const dictionary& dict,
            Random& rndGen
        );


    //- Destructor
    virtual ~fixedValue();


    // Member Functions

         //- Sample the pdf
        virtual scalar sample() const;

        //- Return the minimum value
        virtual scalar minValue() const;

        //- Return the maximum value
        virtual scalar maxValue() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pdfs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
