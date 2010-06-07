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
    Foam::hsPsiThermo

Description
    Sensible enthalpy for a mixture based on compressibility

SourceFiles
    hsPsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hsPsiThermo_H
#define hsPsiThermo_H

#include "basicPsiThermo.H"
#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hsPsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class hsPsiThermo
:
    public basicPsiThermo,
    public MixtureType
{
    // Private data

        //- Sensible enthalpy field [J/kg]
        volScalarField hs_;


    // Private member functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        hsPsiThermo(const hsPsiThermo<MixtureType>&);


public:

    //- Runtime type information
    TypeName("hsPsiThermo");


    // Constructors

        //- Construct from mesh
        hsPsiThermo(const fvMesh&);


    //- Destructor
    virtual ~hsPsiThermo();


    // Member functions

        //- Return the compostion of the mixture
        virtual basicMixture& composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        virtual const basicMixture& composition() const
        {
            return *this;
        }

        //- Update properties
        virtual void correct();


        // Access to thermodynamic state variables

            //- Sensible enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hs()
            {
                return hs_;
            }

            //- Sensible enthalpy [J/kg]
            virtual const volScalarField& hs() const
            {
                return hs_;
            }


        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
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

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;


        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "hsPsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
