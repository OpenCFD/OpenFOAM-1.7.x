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
    Foam::normal

Description
    A normal pdf

    @verbatim
        pdf = strength * exp(-0.5*((x - expectation)/variance)^2 )
    @endverbatim

    strength only has meaning if there's more than one pdf

SourceFiles
    normalI.H
    normal.C
    normalIO.C

\*---------------------------------------------------------------------------*/

#ifndef normal_H
#define normal_H

#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class normal Declaration
\*---------------------------------------------------------------------------*/

class normal
:
    public pdf
{
    // Private data

        dictionary pdfDict_;

        //- min and max values of the distribution
        scalar minValue_;
        scalar maxValue_;

        List<scalar> expectation_;
        List<scalar> variance_;
        List<scalar> strength_;

        scalar range_;


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
