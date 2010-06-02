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

Namespace
    Foam::incompressible::LESModels

Description
    Namespace for incompressible LES models.

Class
    Foam::incompressible::LESModel

Description
    Base class for all incompressible flow LES SGS models.

    This class defines the basic interface for an incompressible flow SGS
    model, and encapsulates data of value to all possible models.
    In particular this includes references to all the dependent fields
    (U, phi), the physical viscosity nu, and the LESProperties
    dictionary, which contains the model selection and model coefficients.

SourceFiles
    LESModel.C

\*---------------------------------------------------------------------------*/

#ifndef LESModel_H
#define LESModel_H

#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "LESdelta.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "incompressible/transportModel/transportModel.H"
#include "wallFvPatch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                           Class LESModel Declaration
\*---------------------------------------------------------------------------*/

class LESModel
:
    public turbulenceModel,
    public IOdictionary
{

protected:

    // Protected data

        Switch printCoeffs_;
        dictionary coeffDict_;

        dimensionedScalar k0_;

        autoPtr<LESdelta> delta_;


    // Protected member functions

        //- Print model coefficients
        virtual void printCoeffs();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        LESModel(const LESModel&);

        //- Disallow default bitwise assignment
        LESModel& operator=(const LESModel&);


public:

    //- Runtime type information
    TypeName("LESModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESModel,
            dictionary,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                transportModel& lamTransportModel
            ),
            (U, phi, lamTransportModel)
        );


    // Constructors

        //- Construct from components
        LESModel
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& lamTransportModel
        );


    // Selectors

        //- Return a reference to the selected LES model
        static autoPtr<LESModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& lamTransportModel
        );


    //- Destructor
    virtual ~LESModel()
    {}


    // Member Functions

        //- Const access to the coefficients dictionary,
        //  which provides info. about choice of models,
        //  and all related data (particularly model coefficients).
        inline const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Access function to filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }

        //- Return the value of k0 which k is not allowed to be less than
        const dimensionedScalar& k0() const
        {
            return k0_;
        }

        //- Allow k0 to be changed
        dimensionedScalar& k0()
        {
            return k0_;
        }


        //- Return the SGS turbulent kinetic energy.
        virtual tmp<volScalarField> k() const = 0;

        //- Return the SGS turbulent dissipation.
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the SGS viscosity.
        virtual tmp<volScalarField> nuSgs() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("nuEff", nuSgs() + nu())
            );
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const = 0;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devBeff() const = 0;

        //- Returns div(dev(Beff)).
        //  This is the additional term due to the filtering of the NSE.
        virtual tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const = 0;


        // RAS compatibility functions for the turbulenceModel base class

            //- Return the turbulence viscosity
            virtual tmp<volScalarField> nut() const
            {
                return nuSgs();
            }

            //- Return the Reynolds stress tensor
            virtual tmp<volSymmTensorField> R() const
            {
                return B();
            }

            //- Return the effective stress tensor including the laminar stress
            virtual tmp<volSymmTensorField> devReff() const
            {
                return devBeff();
            }

            //- Return the source term for the momentum equation
            virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const
            {
                return divDevBeff(U);
            }


        //- Correct Eddy-Viscosity and related properties.
        //  This calls correct(const tmp<volTensorField>& gradU) by supplying
        //  gradU calculated locally.
        void correct();

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
