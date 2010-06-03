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
    Foam::injectorType

Description
    Base class for injectors

SourceFiles
    injectorTypeI.H
    injectorType.C
    injectorTypeIO.C

\*---------------------------------------------------------------------------*/

#ifndef injectorType_H
#define injectorType_H

#include "IOdictionary.H"
#include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "scalarField.H"
#include "Random.H"
#include "liquidMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class injectorType Declaration
\*---------------------------------------------------------------------------*/

class injectorType
{
        typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

public:

    //- Runtime type information
    TypeName("injectorType");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            injectorType,
            dictionary,
            (
                const Time& t,
                const dictionary& dict
            ),
            (t, dict)
        );

    
    // Constructors

        //- Construct from components
        injectorType
        (
            const Time& t,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<injectorType> New
        (
            const Time& t,
            const dictionary& dict
        );


    // Destructor

        virtual ~injectorType();


    // Member Functions
    
        //- Return number of particles to inject
        virtual label nParcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const = 0;

        //- Return the injection position
        virtual const vector position(const label n) const = 0;

        //- Return the injection position
        virtual vector position
        (
            const label n,
            const scalar time,
            const bool twoD,
            const scalar angleOfWedge,
            const vector& axisOfSymmetry,
            const vector& axisOfWedge,
            const vector& axisOfWedgeNormal,
            Random& rndGen
        ) const = 0;

        //- Return the number of holes
        virtual label nHoles() const = 0;
    
        //- Return the injector diameter
        virtual scalar d() const = 0;

        //- Return the injection direction for hole i
        virtual const vector& direction
        (
            const label i, 
            const scalar time
        ) const = 0;

        //- Return the mass of the injected liquid between times
        virtual scalar mass
        (
            const scalar time0,
            const scalar time1,
            const bool twoD,
            const scalar angleOfWedge
        ) const = 0;

        //- Return the mass injected by the injector
        virtual scalar mass() const = 0;

        //- Return the mass flow rate profile
        virtual List<pair> massFlowRateProfile() const = 0;

        //- Return the instantaneous mass flow rate profile
        virtual scalar massFlowRate(const scalar time) const = 0;

        //- Return the pressure injection profile
        virtual List<pair> injectionPressureProfile() const = 0;

        //- Return the instantaneous injection pressure
        virtual scalar injectionPressure(const scalar time) const = 0;

        //- Return the velocity injection profile
        virtual List<pair> velocityProfile() const = 0;

        //- Return the instantaneous velocity
        virtual scalar velocity(const scalar time) const = 0;

        //- Return the discharge coefficient        
        virtual List<pair> CdProfile() const = 0;

        //- Return the instantaneous discharge coefficient
        virtual scalar Cd(const scalar time) const = 0;

        //- Return the fuel mass fractions of the injected particle
        virtual const scalarField& X() const = 0;

        //- Return the temperatue profile of the injected parcel
        virtual List<pair> T() const = 0;

        //- Return the temperatue of the injected parcel
        virtual scalar T(const scalar time) const = 0;

        //- Return the start-of-injection time
        virtual scalar tsoi() const = 0;

        //- Return the end-of-injection time
        virtual scalar teoi() const = 0;

        virtual scalar injectedMass(const scalar t) const = 0;

        virtual bool pressureIndependentVelocity() const = 0;

        //- Return a vector perpendicular to the injection direction and tan2 for hole n
        virtual vector tan1(const label n) const = 0;

        //- Return a vector perpendicular to the injection direction and tan1 for hole n
        virtual vector tan2(const label n) const = 0;

        scalar getTableValue
        (
            const List<pair>& table,
            const scalar value
        ) const;

        scalar integrateTable
        (
            const List<pair>& table,
            const scalar value
        ) const;

        scalar integrateTable
        (
            const List<pair>& table
        ) const;

        virtual void correctProfiles
        (
            const liquidMixture& fuel,
            const scalar referencePressure
        ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
