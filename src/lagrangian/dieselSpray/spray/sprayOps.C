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

#include "spray.H"
#include "atomizationModel.H"
#include "breakupModel.H"
#include "collisionModel.H"
#include "dispersionModel.H"
#include "interpolationCellPoint.H"
#include "processorPolyPatch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void spray::evolve()
{
    sms_.setSize(rho_.size());
    shs_.setSize(rho_.size());
    forAll(srhos_, i)
    {
        srhos_[i].setSize(rho_.size());
    }

    UInterpolator_ = interpolation<vector>::New(interpolationSchemes_, U_);

    rhoInterpolator_ = interpolation<scalar>::New(interpolationSchemes_, rho_);

    pInterpolator_ = interpolation<scalar>::New(interpolationSchemes_, p_);

    TInterpolator_ = interpolation<scalar>::New(interpolationSchemes_, T_);

    calculateAmbientPressure();
    calculateAmbientTemperature();
    collisions().collideParcels(runTime_.deltaT().value());
    move();
    dispersion().disperseParcels();
    inject();
    atomizationLoop();
    breakupLoop();

    UInterpolator_.clear();
    rhoInterpolator_.clear();
    pInterpolator_.clear();
    TInterpolator_.clear();
}


void spray::move()
{
    // Reset Spray Source Terms
    sms_ = vector::zero;
    shs_ = 0.0;
    forAll(srhos_, i)
    {
        srhos_[i] = 0.0;
    }

    Cloud<parcel>::move(*this);
}


void spray::breakupLoop()
{
    forAllIter(spray, *this, elmnt)
    {
        // interpolate...
        vector velocity = UInterpolator().interpolate
        (
            elmnt().position(),
            elmnt().cell()
        );

        // liquidCore < 0.5 indicates discrete drops
        if (elmnt().liquidCore() <= 0.5)
        {
            breakup().updateParcelProperties
            (
                elmnt(),
                runTime_.deltaT().value(),
                velocity,
                fuels_
            );

            breakup().breakupParcel
            (
                elmnt(),
                runTime_.deltaT().value(),
                velocity,
                fuels_
            );
        }
    }
}


void spray::atomizationLoop()
{
    forAllIter(spray, *this, elmnt)
    {
        // interpolate...
        vector velocity = UInterpolator().interpolate
        (
            elmnt().position(),
            elmnt().cell()
        );

        // liquidCore > 0.5 indicates a liquid core
        if (elmnt().liquidCore() > 0.5)
        {
            atomization().atomizeParcel
            (
                elmnt(),
                runTime_.deltaT().value(),
                velocity,
                fuels_
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
