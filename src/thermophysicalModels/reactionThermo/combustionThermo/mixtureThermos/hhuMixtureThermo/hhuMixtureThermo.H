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
    Foam::hhuMixtureThermo

Description
    Foam::hhuMixtureThermo

SourceFiles
    hhuMixtureThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hhuMixtureThermo_H
#define hhuMixtureThermo_H

//#include "hPsiMixtureThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hhuMixtureThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class hhuMixtureThermo
:
    public hhuCombustionThermo,
    public MixtureType
{
    // Private member functions

        void calculate();

        //- Construct as copy (not implemented)
        hhuMixtureThermo(const hhuMixtureThermo<MixtureType>&);


public:

    //- Runtime type information
    TypeName("hhuMixtureThermo");


    // Constructors

        //- Construct from mesh
        hhuMixtureThermo(const fvMesh&);


    //- Destructor
    virtual ~hhuMixtureThermo();


    // Member functions

        //- Return the compostion of the multi-component mixture
        virtual basicMultiComponentMixture& composition()
        {
            return *this;
        }

        //- Return the compostion of the multi-component mixture
        virtual const basicMultiComponentMixture& composition() const
        {
            return *this;
        }

        //- Update properties
        virtual void correct();

        //- Sensible enthalpy [J/kg]
        virtual tmp<volScalarField> hs() const;

        //- Chemical enthalpy [J/kg]
        virtual tmp<volScalarField> hc() const;


        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Unburnt gas enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hu
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Unburnt gas enthalpy for patch [J/kg]
            virtual tmp<scalarField> hu
            (
                const scalarField& T,
                const label patchi
            ) const;


            //- Burnt gas temperature [K]
            virtual tmp<volScalarField> Tb() const;

            //- Unburnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psiu() const;

            //- Burnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psib() const;


        // Access to transport variables

            //- Dynamic viscosity of unburnt gas [kg/ms]
            virtual tmp<volScalarField> muu() const;

            //- Dynamic viscosity of burnt gas [kg/ms]
            virtual tmp<volScalarField> mub() const;


        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "hhuMixtureThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
