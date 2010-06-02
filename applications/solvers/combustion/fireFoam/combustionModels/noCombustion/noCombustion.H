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
    Foam::combustionModel::noCombustion

Description
    No combustion

SourceFiles
    noCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef noCombustion_H
#define noCombustion_H

#include "combustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                          Class noCombustion Declaration
\*---------------------------------------------------------------------------*/

class noCombustion
:
    public combustionModel
{
    // Private data


    // Private Member Functions

        //- Disallow copy construct
        noCombustion(const noCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const noCombustion&);


public:

    //- Runtime type information
    TypeName("noCombustion");


    // Constructors

        //- Construct from components
        noCombustion
        (
            const dictionary& combustionProperties,
            const hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );


    // Destructor

        virtual ~noCombustion();

    // Member Functions

        //- Update properties from given dictionary
        virtual bool read(const dictionary& combustionProperties);

        //- Correct combustion rate
        virtual void correct();

        //- Return normalised consumption rate of (fu - fres)
        virtual tmp<volScalarField> wFuelNorm() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
