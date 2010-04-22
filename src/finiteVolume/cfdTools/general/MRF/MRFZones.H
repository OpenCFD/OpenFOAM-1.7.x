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
    Foam::MRFZones

Description
    Container class for a set of MRFZones with the MRFZone member functions
    implemented to loop over the functions for each MRFZone.

SourceFiles
    MRFZones.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZones_H
#define MRFZones_H

#include "MRFZone.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MRFZones Declaration
\*---------------------------------------------------------------------------*/

class MRFZones
:
    public IOPtrList<MRFZone>
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        MRFZones(const MRFZones&);

        //- Disallow default bitwise assignment
        void operator=(const MRFZones&);


public:

    // Constructors

        //- Construct from fvMesh
        MRFZones(const fvMesh& mesh);


    // Member Functions

        //- Add the Coriolis force contribution to the momentum equation
        void addCoriolis(fvVectorMatrix& UEqn) const;

        //- Add the Coriolis force contribution to the momentum equation
        void addCoriolis(const volScalarField& rho, fvVectorMatrix& UEqn) const;

        //- Make the given absolute velocity relative within the MRF region
        void relativeVelocity(volVectorField& U) const;

        //- Make the given relative velocity absolute within the MRF region
        void absoluteVelocity(volVectorField& U) const;

        //- Make the given absolute flux relative within the MRF region
        void relativeFlux(surfaceScalarField& phi) const;

        //- Make the given absolute mass-flux relative within the MRF region
        void relativeFlux
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Make the given relative flux absolute within the MRF region
        void absoluteFlux(surfaceScalarField& phi) const;

        //- Make the given relative mass-flux absolute within the MRF region
        void absoluteFlux
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Correct the boundary velocity for the roation of the MRF region
        void correctBoundaryVelocity(volVectorField& U) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
