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
    Foam::normal

Description
    A normal pdf

    @verbatim
        pdf = strength * exp(-0.5*((x - expectation)/variance)^2 )
    @endverbatim

    strength only has meaning if there's more than one pdf

SourceFiles
    normal.C

\*---------------------------------------------------------------------------*/

#ifndef normal_H
#define normal_H

#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pdfs
{

/*---------------------------------------------------------------------------*\
                           Class normal Declaration
\*---------------------------------------------------------------------------*/

class normal
:
    public pdf
{
    // Private data


        //- Distribution minimum
        scalar minValue_;

        //- Distribution maximum
        scalar maxValue_;


        // Model coefficients

            scalar expectation_;
            scalar variance_;

            scalar a_;


public:

    //- Runtime type information
    TypeName("normal");


    // Constructors

        //- Construct from components
        normal
        (
            const dictionary& dict,
            Random& rndGen
        );


    //- Destructor
    virtual ~normal();


    // Member Functions

        //- Sample the pdf
        virtual scalar sample() const;

        //- Return the minimum value
        virtual scalar minValue() const;

        //- Return the maximum value
        virtual scalar maxValue() const;

        scalar erfInv(const scalar y) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pdfs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
