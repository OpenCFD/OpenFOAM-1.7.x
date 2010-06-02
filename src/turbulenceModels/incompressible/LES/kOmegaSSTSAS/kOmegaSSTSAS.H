/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::LESmodels::kOmegaSSTSAS

Description
    kOmegaSSTSAS LES turbulence model for incompressible flows

    Reference:
    DESider A European Effort on Hybrid RANS-LES Modelling:
    Results of the European-Union Funded Project, 2004 - 2007
    (Notes on Numerical Fluid Mechanics and Multidisciplinary Design).
    Chapter 8 Formulation of the Scale-Adaptive Simulation (SAS) Model during
    the DESIDER Project.
    F. R. Menter and Y. Egorov.

SourceFiles
    kOmegaSSTSAS.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTSAS_H
#define kOmegaSSTSAS_H

#include "LESModel.H"
#include "volFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class kOmegaSSTSAS Declaration
\*---------------------------------------------------------------------------*/

class kOmegaSSTSAS
:
    public LESModel
{
    // Private member functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volScalarField& D);

        // Disallow default bitwise copy construct and assignment
        kOmegaSSTSAS(const kOmegaSSTSAS&);
        kOmegaSSTSAS& operator=(const kOmegaSSTSAS&);


protected:

    // Protected data

        // Model constants

            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar c1_;
            dimensionedScalar Cs_;

            dimensionedScalar alphaPhi_;
            dimensionedScalar zetaTilda2_;
            dimensionedScalar FSAS_;

            dimensionedScalar omega0_;
            dimensionedScalar omegaSmall_;

            wallDist y_;
            dimensionedScalar Cmu_;
            dimensionedScalar kappa_;


        // Fields

            volScalarField k_;
            volScalarField omega_;
            volScalarField nuSgs_;


    // Protected member functions

        tmp<volScalarField> Lvk2
        (
            const volScalarField& S2
        ) const;

        tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        tmp<volScalarField> F2() const;

        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaK
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        tmp<volScalarField> alphaOmega
        (
             const volScalarField& F1
        ) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        tmp<volScalarField> beta
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, beta1_, beta2_);
        }

        tmp<volScalarField> gamma
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, gamma1_, gamma2_);
        }


public:

    //- Runtime type information
    TypeName("kOmegaSSTSAS");


    // Constructors

        //- Construct from components
        kOmegaSSTSAS
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTSAS()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return omega
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                 new volScalarField("DkEff", alphaK(F1)*nuSgs_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaOmega(F1)*nuSgs_ + nu())
            );
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return SGS viscosity
        virtual tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        virtual tmp<volSymmTensorField> devBeff() const;

        //- Return the deviatoric part of the divergence of Beff
        //  i.e. the additional term in the filtered NSE.
        virtual tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const;

        //- Solve the turbulence equations (k-w) and correct the turbulence
        //  viscosity
        virtual void correct(const tmp<volTensorField>& gradU);

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
