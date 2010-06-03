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
    Foam::blobsSwirlInjector

Description
    Model for parcels injection in pressure Swirl Atomizers.

    Accurate description in:
    @verbatim
        Z. Han, S. Parrish, P.V. Farrell, R.D. Reitz
        "Modeling Atomization Processes Of Pressure Swirl Hollow-Cone Fuel Sprays"
        Atomization and Sprays, vol. 7, pp. 663-684, 1997
    @endverbatim
    and
    @verbatim
        L. Allocca, G. Bella, A. De Vita, L. Di Angelo
        "Experimental Validation of a GDI Spray Model"
        SAE Technical Paper Series, 2002-01-1137
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef blobsSwirlInjector_H
#define blobsSwirlInjector_H

#include "injectorModel.H"
#include "scalarList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class blobsSwirlInjector Declaration
\*---------------------------------------------------------------------------*/

class blobsSwirlInjector
:
    public injectorModel
{

private:

    // Private data

        dictionary blobsSwirlInjectorDict_;

        scalarList coneAngle_;
        scalarList coneInterval_;

        scalarList cD_;
        scalarList cTau_;
        scalarList A_;

    // The initial spray angle for the parcels
        mutable scalar angle_;

    // The initial velocity for the parcels
        mutable scalar u_;

    // The orifice-air cone area ratio

        mutable scalar x_;

    // The initial sheet tickness

        mutable scalar h_;

    // private member functions

        scalar kv
        (
            const label inj
        ) const;

        scalar deltaPressureInj(const scalar time, const label inj) const;
        scalar averagePressure(const label inj) const;

        void calculateHX
        (
            const label inj,
            const scalar massFlow,
            const scalar dPressure,
            const scalar time
        ) const;

public:

    //- Runtime type information
        TypeName("blobsSwirlInjector");


    // Constructors

        //- Construct from components
        blobsSwirlInjector
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~blobsSwirlInjector();


    // Member Functions

        //- Return the injected droplet diameter
        scalar d0(const label injector, const scalar time) const;

        //- Return the spray angle of the injector
        vector direction
        (
            const label injector,
            const label hole,
            const scalar time,
            const scalar d
        ) const;

        scalar velocity
        (
            const label i,
            const scalar time
        ) const;

        scalar averageVelocity
        (
            const label i
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
