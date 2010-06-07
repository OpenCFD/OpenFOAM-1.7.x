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
    Foam::hsReactionThermo

Description
    Sensible enthalpy variant for reactionThermo

SourceFiles
    hReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hsReactionThermo_H
#define hsReactionThermo_H

#include "basicRhoThermo.H"
#include "basicMultiComponentMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class hsReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class hsReactionThermo
:
    public basicRhoThermo
{

protected:

    // Protected data

        //- Sensible enthalpy [J/kg]
        volScalarField hs_;


public:

    //- Runtime type information
    TypeName("hsReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        hsReactionThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from dictionary and mesh
        hsReactionThermo(const fvMesh&);


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<hsReactionThermo> New(const fvMesh&);

        //- Select and check that package contains 'thermoType'
        static autoPtr<hsReactionThermo> NewType
        (
            const fvMesh&,
            const word& thermoType
        );


    //- Destructor
    virtual ~hsReactionThermo();


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


        //- Update properties
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
