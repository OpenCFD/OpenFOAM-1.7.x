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

#include "gradientDispersionRAS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientDispersionRAS, 0);

addToRunTimeSelectionTable
(
    dispersionModel,
    gradientDispersionRAS,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
gradientDispersionRAS::gradientDispersionRAS
(
    const dictionary& dict,
    spray& sm
)
:
    dispersionRASModel(dict, sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gradientDispersionRAS::~gradientDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradientDispersionRAS::disperseParcels() const
{

    const scalar cps = 0.16432;

    scalar dt = spray_.runTime().deltaT().value();
    const volScalarField& k = turbulence().k();
    volVectorField gradk = fvc::grad(k);
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
                vector dir = -gradk[celli]/(mag(gradk[celli]) + SMALL);

                // numerical recipes... Ch. 7. Random Numbers...
                scalar x1 = 0.0;
                scalar x2 = 0.0;
                scalar rsq = 10.0;
                while((rsq > 1.0) || (rsq == 0.0))
                {
                    x1 = 2.0*spray_.rndGen().scalar01() - 1.0;
                    x2 = 2.0*spray_.rndGen().scalar01() - 1.0;
                    rsq = x1*x1 + x2*x2;
                }
                
                scalar fac = sqrt(-2.0*log(rsq)/rsq);
                
                // in 2D calculations the -grad(k) is always
                // away from the axis of symmetry
                // This creates a 'hole' in the spray and to
                // prevent this we let x1 be both negative/positive
                if (spray_.twoD())
                {
                    fac *= x1;
                }
                else
                {
                    fac *= mag(x1);
                }
                
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
