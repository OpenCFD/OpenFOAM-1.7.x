/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::speciesTransport

Description
    Basic species transport type based on the use of a fitting
    function for nu.

    All other properties are derived from this primitive function.

    The nu function must be provided by the derived type,
    e.g. SutherlandTransport.

SourceFiles
    speciesTransportI.H
    speciesTransport.C


\*---------------------------------------------------------------------------*/

#ifndef speciesTransport_H
#define speciesTransport_H

//#include "speciesThermo.H"
#include "janafThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class speciesTransport Declaration
\*---------------------------------------------------------------------------*/

class speciesTransport
:
    public janafThermo
{

public:

    // Constructors

        //- Construct from speciesThermo
        inline speciesTransport
        (
            const janafThermo& sThermo
        );

        //- Construct from Istream
        speciesTransport(Istream&);


    // Member Functions

        // Dynamic viscosity [kg/ms]
        //inline scalar mu(const scalar T) const

        // Thermal conductivity [W/mK]
        //inline scalar kappa(const scalar T) const;

        // Thermal diffusivity for enthalpy [kg/ms]
        //inline scalar alpha(const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar T) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const speciesTransport&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "speciesTransportI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
