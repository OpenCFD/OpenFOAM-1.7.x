/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::particleForces

Description
    Provides a mechanism to calculate particle forces
    Note: forces are force per unit mass (accelerations)

SourceFiles
    particleForces.C

\*---------------------------------------------------------------------------*/

#ifndef particleForces_H
#define particleForces_H

#include "dictionary.H"
#include "Switch.H"
#include "vector.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class particleForces Declaration
\*---------------------------------------------------------------------------*/

class particleForces
{
    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- The particleForces dictionary
        const dictionary& dict_;

        //- Gravity
        const vector g_;

        //- Velocity gradient field
        const volTensorField* gradUPtr_;


        // Forces to include in particle motion evaluation

            //- Gravity
            Switch gravity_;

            //- Virtual mass
            Switch virtualMass_;

            //- Virtual mass force coefficient
            scalar Cvm_;

            //- Pressure gradient
            Switch pressureGradient_;


        // Additional info

            //- Name of velucity field - default = "U"
            const word UName_;


public:

    // Constructors

        //- Construct from mesh, dictionary and gravity
        particleForces
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const vector& g
        );

        //- Construct copy
        particleForces(const particleForces& f);


    //- Destructor
    ~particleForces();


    // Member Functions

        // Access

            //- Return the particleForces dictionary
            const dictionary& dict() const;

            //- Return the gravity vector
            const vector& g() const;

            //- Return gravity force activate switch
            Switch gravity() const;

            //- Return virtual mass force activate switch
            Switch virtualMass() const;

            //- Return virtual mass force coefficient
            Switch Cvm() const;

            //- Return pressure gradient force activate switch
            Switch pressureGradient() const;

            //- Return name of velocity field
            const word& UName() const;


       // Evaluation

            //- Cache carrier fields
            void cacheFields(const bool store);

            //- Calculate action/reaction forces between carrier and particles
            vector calcCoupled
            (
                const label cellI,
                const scalar dt,
                const scalar rhoc,
                const scalar rho,
                const vector& Uc,
                const vector& U
            ) const;

            //- Calculate external forces applied to the particles
            vector calcNonCoupled
            (
                const label cellI,
                const scalar dt,
                const scalar rhoc,
                const scalar rho,
                const vector& Uc,
                const vector& U
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
