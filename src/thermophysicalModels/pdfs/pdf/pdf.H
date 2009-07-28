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
    Foam::pdf

Description
    A library of runtime-selectable pdf's.

    Returns a sampled value given the expectation (nu) and variance (sigma^2)

    Sample of planned pdfs (beta p. 374-375):
    - binomial
    - geometric
    - Poisson
    - hypergeometric
    - Pascal
    - uniform
    - exponential
    - normal
    - Gamma
    - chi
    - t?
    - F?
    - beta
    - Weibull
    - Rayleigh
    - Cauchy?

    The pdf is tabulated in equidistant nPoints, in an interval.
    These values are integrated to obtain the cumulated pdf,
    which is then used to change the distribution from unifrom to
    the actual pdf.

SourceFiles
    pdfI.H
    pdf.C
    pdfIO.C

\*---------------------------------------------------------------------------*/

#ifndef pdf_H
#define pdf_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdf Declaration
\*---------------------------------------------------------------------------*/

class pdf
{

protected:

    // Protected data

        const dictionary& dict_;

        Random& rndGen_;


public:

    //-Runtime type information
    TypeName("pdf");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        pdf,
        dictionary,
        (
            const dictionary& dict,
            Random& rndGen
        ),
        (dict, rndGen)
    );


    // Constructors

        //- Construct from dictionary
        pdf(const dictionary& dict, Random& rndGen);


    //- Selector
    static autoPtr<pdf> New(const dictionary& dict, Random& rndGen);


    //- Destructor
    virtual ~pdf();


    // Member Functions

        //- Sample the pdf
        virtual scalar sample() const = 0;

        //- Return the minimum value
        virtual scalar minValue() const = 0;

        //- Return the maximum value
        virtual scalar maxValue() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
