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
    Foam::compressible::RASModels::PDRkEpsilon

Description
    Standard k-epsilon turbulence model with additional source terms
    corresponding to PDR basic drag model (\link basic.H \endlink)

    The default model coefficients correspond to the following:
    @verbatim
        kEpsilonCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          -0.33;  // only for compressible
            sigmak      1.0;    // only for compressible
            sigmaEps    1.3;
            Prt         1.0;    // only for compressible
        }
    @endverbatim

    The turbulence source term \f$ G_{R} \f$ appears in the
    \f$ \kappa-\epsilon \f$ equation for the generation of turbulence due to
    interaction with unresolved obstacles.

    In the \f$ \epsilon  \f$ equation \f$ C_{1} G_{R} \f$ is added as a source
    term.

    In the \f$ \kappa \f$ equation \f$ G_{R} \f$ is added as a source term.

SourceFiles
    PDRkEpsilon.C
    PDRkEpsilonCorrect.C

\*---------------------------------------------------------------------------*/

#ifndef compressiblePDRkEpsilon_H
#define compressiblePDRkEpsilon_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class PDRkEpsilon Declaration
\*---------------------------------------------------------------------------*/

class PDRkEpsilon
:
    public RASModel
{
    // Private data

    // Model coefficients

        dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar sigmak_;
        dimensionedScalar sigmaEps_;
        dimensionedScalar Prt_;

    // Fields

        volScalarField k_;
        volScalarField epsilon_;
        volScalarField mut_;
        volScalarField alphat_;


public:

    //- Runtime type information
    TypeName("PDRkEpsilon");


    // Constructors

        //- Construct from components
        PDRkEpsilon
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermophysicalModel
        );


    //- Destructor
    virtual ~PDRkEpsilon()
    {}


    // Member Functions

        tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", mut_/sigmak_ + mu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", mut_/sigmaEps_ + mu())
            );
        }

        //- Return the effective turbulent thermal diffusivity
        tmp<volScalarField> alphaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("alphaEff", alphat_ + alpha())
            );
        }

        //- Return the turbulence kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
