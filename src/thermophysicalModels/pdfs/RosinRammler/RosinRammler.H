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
    Foam::RosinRammler

Description
    Rosin-Rammler pdf

   @f[
       cumulative pdf = (1.0 - exp( -((x - d0)/d)^n ) / (1.0 - exp( -((d1 - d0)/d)^n )
   @f]


SourceFiles
    RosinRammler.C

\*---------------------------------------------------------------------------*/

#ifndef RosinRammler_H
#define RosinRammler_H

#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pdfs
{

/*---------------------------------------------------------------------------*\
                       Class RosinRammler Declaration
\*---------------------------------------------------------------------------*/

class RosinRammler
:
    public pdf
{
    // Private data

        //- Distribution minimum
        scalar minValue_;

        //- Distribution maximum
        scalar maxValue_;

        // Model coefficients

            scalar d_;
            scalar n_;


public:

    //- Runtime type information
    TypeName("RosinRammler");


    // Constructors

        //- Construct from components
        RosinRammler
        (
            const dictionary& dict,
            Random& rndGen
        );


    //- Destructor
    virtual ~RosinRammler();


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
