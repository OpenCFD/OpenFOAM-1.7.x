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
    Foam::Gibilaro

Description
    H, Enwald, E. Peirano, A-E Almstedt
    'Eulerian Two-Phase Flow Theory Applied to Fluidization'
    Int. J. Multiphase Flow, Vol. 22, Suppl, pp. 21-66 (1996)
    Eq. 106, p. 43

SourceFiles
    Gibilaro.C

\*---------------------------------------------------------------------------*/

#ifndef Gibilaro_H
#define Gibilaro_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Gibilaro Declaration
\*---------------------------------------------------------------------------*/

class Gibilaro
:
    public dragModel
{

public:

    //- Runtime type information
    TypeName("Gibilaro");


    // Constructors

        //- Construct from components
        Gibilaro
        (
            const dictionary& interfaceDict,
            const volScalarField& alpha,
            const phaseModel& phasea,
            const phaseModel& phaseb
        );


    //- Destructor
    virtual ~Gibilaro();


    // Member Functions

        tmp<volScalarField> K(const volScalarField& Ur) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
