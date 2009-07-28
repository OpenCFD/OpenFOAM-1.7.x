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
    Foam::compressibilityModels::Wallis

Description
    Wallis compressibility model.

SourceFiles
    Wallis.C

\*---------------------------------------------------------------------------*/

#ifndef Wallis_H
#define Wallis_H

#include "barotropicCompressibilityModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressibilityModels
{

/*---------------------------------------------------------------------------*\
                           Class Wallis Declaration
\*---------------------------------------------------------------------------*/

class Wallis
:
    public barotropicCompressibilityModel
{
    // Private data

        dimensionedScalar psiv_;
        dimensionedScalar psil_;

        dimensionedScalar rhovSat_;
        dimensionedScalar rholSat_;


public:

    //- Runtime type information
    TypeName("Wallis");


    // Constructors

        //- construct from components
        Wallis
        (
            const dictionary& compressibilityProperties,
            const volScalarField& gamma,
            const word& psiName = "psi"
        );


    // Destructor

        ~Wallis()
        {}


    // Member Functions

        //- Correct the Wallis compressibility
        void correct();

        //- Read transportProperties dictionary
        bool read(const dictionary& compressibilityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressibilityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
