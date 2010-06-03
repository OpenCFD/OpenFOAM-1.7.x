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

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "stochasticDispersionRAS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(stochasticDispersionRAS, 0);

addToRunTimeSelectionTable
(
    dispersionModel,
    stochasticDispersionRAS,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
stochasticDispersionRAS::stochasticDispersionRAS
(
    const dictionary& dict,
    spray& sm
)
:
    dispersionRASModel(dict, sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

stochasticDispersionRAS::~stochasticDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void stochasticDispersionRAS::disperseParcels() const
{

    const scalar cps = 0.16432;
    const vector one(1.0, 1.0, 1.0);

    scalar dt = spray_.runTime().deltaT().value();
    const volScalarField& k = turbulence().k();
    //volVectorField gradk = fvc::grad(k);
    const volScalarField& epsilon = turbulence().epsilon();
    const volVectorField& U = spray_.U();

    for
    (
        spray::iterator elmnt = spray_.begin();
        elmnt != spray_.end();
        ++elmnt
    )
    {
        label celli = elmnt().cell();
        scalar UrelMag = mag(elmnt().U() - U[celli] - elmnt().Uturb());

        scalar Tturb = min
        (
            k[celli]/epsilon[celli], 
            cps*pow(k[celli], 1.5)/epsilon[celli]/(UrelMag + SMALL)
        );

        // parcel is perturbed by the turbulence
        if (dt < Tturb)
        {
            elmnt().tTurb() += dt;

            if (elmnt().tTurb() > Tturb)
            {
                elmnt().tTurb() = 0.0;
                
                scalar sigma = sqrt(2.0*k[celli]/3.0);
                vector dir = 2.0*spray_.rndGen().vector01() - one;
                dir /= mag(dir) + SMALL;
                
                // numerical recipes... Ch. 7. Random Numbers...
                scalar x1,x2;
                scalar rsq = 10.0;
                while((rsq > 1.0) || (rsq == 0.0))
                {
                    x1 = 2.0*spray_.rndGen().scalar01() - 1.0;
                    x2 = 2.0*spray_.rndGen().scalar01() - 1.0;
                    rsq = x1*x1 + x2*x2;
                }
                
                scalar fac = sqrt(-2.0*log(rsq)/rsq);
                
                fac *= mag(x1);
                
                elmnt().Uturb() = sigma*fac*dir;

            }
        }
        else
        {
            elmnt().tTurb() = GREAT;
            elmnt().Uturb() = vector::zero;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
