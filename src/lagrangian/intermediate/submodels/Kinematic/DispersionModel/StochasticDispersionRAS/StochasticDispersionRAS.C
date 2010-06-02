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

\*---------------------------------------------------------------------------*/

#include "StochasticDispersionRAS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticDispersionRAS<CloudType>::StochasticDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticDispersionRAS<CloudType>::~StochasticDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::StochasticDispersionRAS<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::vector Foam::StochasticDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    const scalar cps = 0.16432;

    const volScalarField& k = *this->kPtr_;
    const volScalarField& epsilon = *this->epsilonPtr_;

    const scalar UrelMag = mag(U - Uc - UTurb);

    const scalar tTurbLoc = min
    (
        k[celli]/epsilon[celli],
        cps*pow(k[celli], 1.5)/epsilon[celli]/(UrelMag + SMALL)
    );

    // Parcel is perturbed by the turbulence
    if (dt < tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0.0;

            scalar sigma = sqrt(2.0*k[celli]/3.0);
            vector dir = 2.0*this->owner().rndGen().vector01() - vector::one;
            dir /= mag(dir) + SMALL;

            // Numerical Recipes... Ch. 7. Random Numbers...
            scalar x1 = 0.0;
            scalar x2 = 0.0;
            scalar rsq = 10.0;
            while ((rsq > 1.0) || (rsq == 0.0))
            {
                x1 = 2.0*this->owner().rndGen().scalar01() - 1.0;
                x2 = 2.0*this->owner().rndGen().scalar01() - 1.0;
                rsq = x1*x1 + x2*x2;
            }

            scalar fac = sqrt(-2.0*log(rsq)/rsq);

            fac *= mag(x1);

            UTurb = sigma*fac*dir;

        }
    }
    else
    {
        tTurb = GREAT;
        UTurb = vector::zero;
    }

    return Uc + UTurb;
}


// ************************************************************************* //
