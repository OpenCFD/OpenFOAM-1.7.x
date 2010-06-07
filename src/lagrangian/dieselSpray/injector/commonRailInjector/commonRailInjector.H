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
    Foam::commonRailInjector

Description
    Common-rail injector

SourceFiles
    commonRailInjectorI.H
    commonRailInjector.C
    commonRailInjectorIO.C

\*---------------------------------------------------------------------------*/

#ifndef commonRailInjector_H
#define commonRailInjector_H

#include "injectorType.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class commonRailInjector Declaration
\*---------------------------------------------------------------------------*/

class commonRailInjector
:
    public injectorType
{

private:

    typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

    // Private data

        dictionary propsDict_;

        vector position_;
        vector direction_;
        scalar d_;
        scalar mass_;
        scalar injectionPressure_;
        scalar T_;
        label nParcels_;
        scalarField X_;
        List<pair> massFlowRateProfile_;
        List<pair> velocityProfile_;
        List<pair> injectionPressureProfile_;
        List<pair> CdProfile_;
        List<pair> TProfile_;
        scalar averageParcelMass_;

        bool pressureIndependentVelocity_;

        //- two orthogonal vectors that are also orthogonal
        //  to the injection direction
        vector tangentialInjectionVector1_, tangentialInjectionVector2_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        commonRailInjector(const commonRailInjector&);

        //- Disallow default bitwise assignment
        void operator=(const commonRailInjector&);

        //- Create two vectors orthonoal to each other
        //  and the injection vector
        void setTangentialVectors();

        //- Return the fraction of the total injected liquid
        scalar fractionOfInjection(const scalar time) const;


public:

    //- Runtime type information
    TypeName("commonRailInjector");


    // Constructors

        //- Construct from components
        commonRailInjector
        (
            const Time& t,
            const dictionary& dict
        );


    // Destructor

        ~commonRailInjector();


    // Member Functions

        //- Return number of particles to inject
        label nParcelsToInject
        (
            const scalar t0,
            const scalar t1
        ) const;

        //- Return the injection position
        const vector position(const label n) const;

        //- Return the injection position
        vector position
        (
            const label n,
            const scalar time,
            const bool twoD,
            const scalar angleOfWedge,
            const vector& axisOfSymmetry,
            const vector& axisOfWedge,
            const vector& axisOfWedgeNormal,
            Random& rndGen
        ) const;
    
        //- Return the number of holes
        label nHoles() const;

        //- Return the injector diameter
        scalar d() const;

        //- Return the injection direction
        const vector& direction
        (
            const label i,
            const scalar time
        ) const;

        //- Return the mass of the injected particle
        scalar mass
        (
            const scalar t0,
            const scalar t1,
            const bool twoD,
            const scalar angleOfWedge
        ) const;

        //- Return the mass injected by the injector
        scalar mass() const;

        //- Return the fuel mass fractions of the injected particle
        const scalarField& X() const;

        //- Return the temperature profile of the injected particle
        List<pair> T() const;

        //- Return the temperature of the injected particle
        scalar T(const scalar time) const;

        //- Return the start-of-injection time
        scalar tsoi() const;

        //- Return the end-of-injection time
        scalar teoi() const;

        //- Return the injected liquid mass
        scalar injectedMass(const scalar t) const;

        List<pair> massFlowRateProfile() const
        {
            return massFlowRateProfile_;
        }

        scalar massFlowRate(const scalar time) const;

        List<pair> injectionPressureProfile() const
        {
            return injectionPressureProfile_;
        }
    
        scalar injectionPressure(const scalar time) const;

        List<pair> velocityProfile() const
        {
            return velocityProfile_;
        }
        
        scalar velocity(const scalar time) const;

        List<pair> CdProfile() const;
        scalar Cd(const scalar time) const;

        vector tan1(const label n) const;
        vector tan2(const label n) const;

        void correctProfiles
        (
            const liquidMixture& fuel,
            const scalar referencePressure
        );

        bool pressureIndependentVelocity() const
        {
            return pressureIndependentVelocity_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
