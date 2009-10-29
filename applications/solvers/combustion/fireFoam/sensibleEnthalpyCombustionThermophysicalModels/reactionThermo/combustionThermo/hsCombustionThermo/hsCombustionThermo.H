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
    Foam::hsCombustionThermo

Description
    Foam::hsCombustionThermo

SourceFiles
    hsCombustionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hsCombustionThermo_H
#define hsCombustionThermo_H

#include "basicSensiblePsiThermo.H"
#include "basicMultiComponentMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class hsCombustionThermo Declaration
\*---------------------------------------------------------------------------*/

class hsCombustionThermo
:
    public basicSensiblePsiThermo
{

protected:

    // Protected data

        //- sensible enthalpy field
        volScalarField hs_;


public:

    //- Runtime type information
    TypeName("hsCombustionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        hsCombustionThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from dictionary and mesh
        hsCombustionThermo(const fvMesh&);


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<hsCombustionThermo> New(const fvMesh&);

        //- Select and check that package contains 'thermoType'
        static autoPtr<hsCombustionThermo> NewType
        (
            const fvMesh&,
            const word& thermoType
        );


    //- Destructor
    virtual ~hsCombustionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicMultiComponentMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicMultiComponentMixture& composition() const = 0;


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

        //- Chemical enthalpy [J/kg]
        virtual tmp<volScalarField> hc() const = 0;

        //- Update properties
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
