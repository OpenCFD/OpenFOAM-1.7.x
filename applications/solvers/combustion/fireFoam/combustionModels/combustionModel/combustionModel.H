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
    Foam::combustionModel

Description
    Base class for all non-premixed combustion models.

SourceFiles
    combustionModel.C

\*---------------------------------------------------------------------------*/

#ifndef combustionModel_H
#define combustionModel_H

#include "IOdictionary.H"
#include "hsCombustionThermo.H"
#include "turbulenceModel.H"
#include "multivariateSurfaceInterpolationScheme.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class combustionModel Declaration
\*---------------------------------------------------------------------------*/

class combustionModel
{

protected:

    // Protected data

        //- Dictionary of coefficients for the particular model
        dictionary combustionModelCoeffs_;

        //- Reference to the thermodynamic
        const hsCombustionThermo& thermo_;

        //- Reference to the turbulence model
        const compressible::turbulenceModel& turbulence_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to mass-flux field
        const surfaceScalarField& phi_;

        //- Reference to the density field
        const volScalarField& rho_;

        //- Stoichiometric air-fuel mass ratio
        dimensionedScalar stoicRatio_;

        //- Stoichiometric oxygen-fuel mass ratio
        dimensionedScalar s_;

        //- Heat of combustion (J/Kg)
        dimensionedScalar qFuel_;


private:

    // Private Member Functions

        //- Disallow copy construct
        combustionModel(const combustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const combustionModel&);

        const basicMultiComponentMixture& composition_;


public:

    //- Runtime type information
    TypeName("combustionModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            combustionModel,
            dictionary,
            (
                const dictionary& combustionProperties,
                const hsCombustionThermo& thermo,
                const compressible::turbulenceModel& turbulence,
                const surfaceScalarField& phi,
                const volScalarField& rho
            ),
            (
                combustionProperties,
                thermo,
                turbulence,
                phi,
                rho
            )
        );


    // Selectors

        //- Return a reference to the selected combustion model
        static autoPtr<combustionModel> New
        (
            const dictionary& combustionProperties,
            const hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );


    // Constructors

        //- Construct from components
        combustionModel
        (
            const dictionary& combustionProperties,
            const hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );


    //- Destructor
    virtual ~combustionModel();


    // Member Functions

        // Access functions

            //- Access composition
            const basicMultiComponentMixture& composition() const
            {
                return composition_;
            }

            //- Access combustion dictionary
            const dictionary combustionModelCoeffs() const
            {
                return combustionModelCoeffs_;
            }

            //- Access heat of combustion
            const dimensionedScalar qFuel() const
            {
                return qFuel_;
            }

        //- Return normalised consumption rate of (fu - fres)
        virtual tmp<volScalarField> wFuelNorm() const = 0;

        //- Fuel consumption rate matrix i.e. source-term for the fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& fu) const;

        //- Heat-release rate calculated from the given
        //  fuel consumption rate matrix
        virtual tmp<volScalarField> dQ(const fvScalarMatrix& Rfu) const;

        //- Correct combustion rate
        virtual void correct() = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& combustionProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
