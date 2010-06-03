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
    Foam::incompressible::LESModels::Smagorinsky2

Description
    The Isochoric Smagorinsky Model for incompressible flows

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevails, hence
    @verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D) - 2*cD2*delta*(D.dev(D));
        Beff = 2/3*k*I - 2*nuEff*dev(D) - 2*cD2*delta*(D.dev(D));

    where

        D = symm(grad(U));
        k = cI*delta^2*||D||^2
        nuSgs = ck*sqrt(k)*delta
        nuEff = nuSgs + nu
    @endverbatim

SourceFiles
    Smagorinsky2.C

\*---------------------------------------------------------------------------*/

#ifndef Smagorinsky2_H
#define Smagorinsky2_H

#include "Smagorinsky.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class Smagorinsky2 Declaration
\*---------------------------------------------------------------------------*/

class Smagorinsky2
:
    public Smagorinsky
{
    // Private data

        dimensionedScalar cD2_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        Smagorinsky2(const Smagorinsky2&);
        Smagorinsky2& operator=(const Smagorinsky2&);


public:

    //- Runtime type information
    TypeName("Smagorinsky2");


    // Constructors

        //- Construct from components
        Smagorinsky2
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~Smagorinsky2()
    {}


    // Member Functions

        //- Return B.
        virtual tmp<volSymmTensorField> B() const;

        //- Returns div(B).
        // This is the additional term due to the filtering of the NSE.
        virtual tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const;

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
