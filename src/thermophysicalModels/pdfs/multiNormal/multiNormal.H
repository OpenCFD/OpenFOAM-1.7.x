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
    Foam::multiNormal

Description
    A multiNormal pdf

    @verbatim
        pdf = sum_i strength_i * exp(-0.5*((x - expectation_i)/variance_i)^2 )
    @endverbatim


SourceFiles
    multiNormal.C

\*---------------------------------------------------------------------------*/

#ifndef multiNormal_H
#define multiNormal_H

#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pdfs
{

/*---------------------------------------------------------------------------*\
                        Class multiNormal Declaration
\*---------------------------------------------------------------------------*/

class multiNormal
:
    public pdf
{
    // Private data

        //- Distribution minimum
        scalar minValue_;

        //- Distribution maximum
        scalar maxValue_;

        //- Distribution range
        scalar range_;


        // Model coefficients

            List<scalar> expectation_;
            List<scalar> variance_;
            List<scalar> strength_;


public:

    //- Runtime type information
    TypeName("multiNormal");


    // Constructors

        //- Construct from components
        multiNormal
        (
            const dictionary& dict,
            Random& rndGen
        );


    //- Destructor
    virtual ~multiNormal();


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
