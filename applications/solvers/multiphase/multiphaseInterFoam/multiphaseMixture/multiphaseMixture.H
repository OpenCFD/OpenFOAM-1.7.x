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
    Foam::multiphaseMixture

Description
    Incompressible multi-phase mixture with built in solution for the
    phase fractions with interface compression for interface-capturing.

    Derived from transportModel so that it can be unsed in conjunction with
    the incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface
    between each phase-pair.

SourceFiles
    multiphaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseMixture_H
#define multiphaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "phase.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "multivariateSurfaceInterpolationScheme.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class multiphaseMixture Declaration
\*---------------------------------------------------------------------------*/

class multiphaseMixture
:
    public transportModel
{
public:

    class interfacePair
    :
        public Pair<word>
    {
    public:

        class hash
        :
            public Hash<interfacePair>
        {
        public:

            hash()
            {}

            label operator()(const interfacePair& key) const
            {
                return word::hash()(key.first()) + word::hash()(key.second());
            }
        };


        // Constructors

            interfacePair()
            {}

            interfacePair(const word& alpha1Name, const word& alpha2Name)
            :
                Pair<word>(alpha1Name, alpha2Name)
            {}

            interfacePair(const phase& alpha1, const phase& alpha2)
            :
                Pair<word>(alpha1.name(), alpha2.name())
            {}


        // Friend Operators

            friend bool operator==
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return
                (
                    ((a.first() == b.first()) && (a.second() == b.second()))
                 || ((a.first() == b.second()) && (a.second() == b.first()))
                );
            }

            friend bool operator!=
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (!(a == b));
            }
    };


private:

    // Private data

        //- Dictionary of phases
        PtrDictionary<phase> phases_;

        //- The phase chosen as reference, the one which is derived from
        //  the others such thatr they sum to 1
        phase& refPhase_;

        const fvMesh& mesh_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        surfaceScalarField rhoPhi_;

        volScalarField alphas_;

        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            sigmaTable;

        sigmaTable sigmas_;
        dimensionSet dimSigma_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- Conversion factor for degrees into radians
        static const scalar convertToRad;

        //- Phase-fraction field table for multivariate discretisation
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable alphaTable_;


    // Private member functions

        void calcAlphas();

        void solveAlphas
        (
            const label nAlphaCorr,
            const bool cycleAlpha,
            const scalar cAlpha
        );

        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        void correctContactAngle
        (
            const phase& alpha1,
            const phase& alpha2,
            surfaceVectorField::GeometricBoundaryField& nHatb
        ) const;

        tmp<volScalarField> K(const phase& alpha1, const phase& alpha2) const;


public:

    // Constructors

        //- Construct from components
        multiphaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        ~multiphaseMixture()
        {}


    // Member Functions

        //- Return the phases
        const PtrDictionary<phase>& phases() const
        {
            return phases_;
        }

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        const surfaceScalarField& rhoPhi() const
        {
            return rhoPhi_;
        }

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Correct the mixture properties
        void correct();

        //- Read base transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
