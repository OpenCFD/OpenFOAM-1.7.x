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
    Foam::compressible::LESModels::dynOneEqEddy

Description
    One Equation Eddy Viscosity Model for compressible flows.

    Eddy viscosity SGS model using a modeled balance equation to simulate
    the behaviour of k.
    Thus
    @verbatim
        d/dt(k) + div(U*k) - div(nuSgs*grad(k))
        =
        -rho*B*L - ce*rho*k^3/2/delta

    and

        B = 2/3*k*I - 2*nuSgs*dev(D)

    where

        D = symm(grad(U));
        nuSgs = ck*sqrt(k)*delta
    @endverbatim

SourceFiles
    dynOneEqEddy.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleDynOneEqEddy_H
#define compressibleDynOneEqEddy_H

#include "GenEddyVisc.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class dynOneEqEddy Declaration
\*---------------------------------------------------------------------------*/

class dynOneEqEddy
:
    public GenEddyVisc
{
    // Private data

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volSymmTensorField& D);

        //- Calculate ck, ce by filtering the velocity field U.
        dimensionedScalar ck_(const volSymmTensorField& D) const;
        dimensionedScalar ce_(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        dynOneEqEddy(const dynOneEqEddy&);
        dynOneEqEddy& operator=(const dynOneEqEddy&);


public:

    //- Runtime type information
    TypeName("dynOneEqEddy");


    // Constructors

        //- Constructor from components
        dynOneEqEddy
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
        );


    // Destructor
    virtual ~dynOneEqEddy()
    {}


    // Member Functions

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", muSgs_ + mu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
