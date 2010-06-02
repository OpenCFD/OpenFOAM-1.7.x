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

#include "spray.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar spray::injectedMass(const scalar t) const
{
    scalar sum = 0.0;

    forAll (injectors_, i)
    {
        sum += injectors_[i].properties()->injectedMass(t);
    }

    return sum;
}


scalar spray::totalMassToInject() const
{
    scalar sum = 0.0;

    forAll (injectors_, i)
    {
        sum += injectors_[i].properties()->mass();
    }

    return sum;
}


scalar spray::injectedEnthalpy
(
    const scalar time
) const
{
    scalar sum = 0.0;
    label Nf = fuels_->components().size();

    forAll (injectors_, i)
    {
        scalar T = injectors_[i].properties()->T(time);
        scalarField X(injectors_[i].properties()->X());
        scalar pi = 1.0e+5;
        scalar hl = fuels_->hl(pi, T, X);
        scalar Wl = fuels_->W(X);
        scalar hg = 0.0;

        for(label j=0; j<Nf; j++)
        {
            label k = liquidToGasIndex_[j];
            hg += gasProperties()[k].H(T)*gasProperties()[k].W()*X[j]/Wl;
        }

        sum += injectors_[i].properties()->injectedMass(time)*(hg-hl);
    }

    return sum;
}


scalar spray::liquidMass() const
{
    scalar sum = 0.0;

    for
    (
        spray::const_iterator elmnt = begin();
        elmnt != end();
        ++elmnt
    )
    {
        sum += elmnt().m();
    }

    if (twoD())
    {
        sum *= 2.0*mathematicalConstant::pi/angleOfWedge();
    }

    reduce(sum, sumOp<scalar>());

    return sum;
}


scalar spray::liquidEnthalpy() const
{
    scalar sum = 0.0;
    label Nf = fuels().components().size();

    for
    (
        spray::const_iterator elmnt = begin();
        elmnt != end();
        ++elmnt
    )
    {
        scalar T = elmnt().T();
        scalar pc = p()[elmnt().cell()];
        scalar hlat = fuels().hl(pc, T, elmnt().X());
        scalar hg = 0.0;
        scalar Wl = fuels().W(elmnt().X());

        for(label j=0; j<Nf; j++)
        {
            label k = liquidToGasIndex_[j];

            hg += 
                gasProperties()[k].H(T)*gasProperties()[k].W()*elmnt().X()[j]
               /Wl;
        }

        scalar h = hg - hlat;
        sum += elmnt().m()*h;
    }

    if (twoD())
    {
        sum *= 2.0*mathematicalConstant::pi/angleOfWedge();
    }

    reduce(sum, sumOp<scalar>());

    return sum;
}


scalar spray::liquidTotalEnthalpy() const
{
    scalar sum = 0.0;
    label Nf = fuels().components().size();

    for
    (
        spray::const_iterator elmnt = begin();
        elmnt != end();
        ++elmnt
    )
    {
        label celli = elmnt().cell();
        scalar T = elmnt().T();
        scalar pc = p()[celli];
        scalar rho = fuels().rho(pc, T, elmnt().X());
        scalar hlat = fuels().hl(pc, T, elmnt().X());
        scalar hg = 0.0;
        scalar Wl = fuels().W(elmnt().X());

        for(label j=0; j<Nf; j++)
        {
            label k = liquidToGasIndex_[j];
            hg += 
                gasProperties()[k].H(T)*gasProperties()[k].W()*elmnt().X()[j]
               /Wl;
        }

        scalar psat = fuels().pv(pc, T, elmnt().X());

        scalar h = hg - hlat + (pc - psat)/rho;
        sum += elmnt().m()*h;
    }

    if (twoD())
    {
        sum *= 2.0*mathematicalConstant::pi/angleOfWedge();
    }

    reduce(sum, sumOp<scalar>());

    return sum;
}


scalar spray::liquidKineticEnergy() const
{
    scalar sum = 0.0;
    for
    (
        spray::const_iterator elmnt = begin();
        elmnt != end();
        ++elmnt
    )
    {
        scalar ke = pow(mag(elmnt().U()), 2.0);
        sum += elmnt().m()*ke;
    }

    if (twoD())
    {
        sum *= 2.0*mathematicalConstant::pi/angleOfWedge();
    }

    reduce(sum, sumOp<scalar>());

    return 0.5*sum;

}


scalar spray::injectedLiquidKineticEnergy() const
{
    return injectedLiquidKE_;
}


scalar spray::liquidPenetration(const scalar prc) const
{
    return liquidPenetration(0, prc);
}


scalar spray::liquidPenetration
(
    const label nozzlei,
    const scalar prc
) const
{

    label nHoles = injectors_[nozzlei].properties()->nHoles();
    vector ip(vector::zero);
    if (nHoles > 1)
    {
        for(label i=0;i<nHoles;i++)
        {
            ip += injectors_[nozzlei].properties()->position(i);
        }
        ip /= nHoles;
    }
    else
    {
        ip = injectors_[nozzlei].properties()->position(0);
    }

//    vector ip = injectors_[nozzlei].properties()->position();
    scalar d = 0.0;
    scalar mTot = 0.0;

    label Np = size();
    
    // arrays containing the parcels mass and
    // distance from injector in ascending order
    scalarField m(Np);
    scalarField dist(Np);
    label n = 0;

    if (Np > 1)
    {
        // NN.
        // first arrange the parcels in ascending order
        // the first parcel is closest to injector
        // and the last one is most far away.
        spray::const_iterator first = begin();
        m[n] = first().m();
        dist[n] = mag(first().position() - ip);

        mTot += m[n];

        for
        (
            spray::const_iterator elmnt = ++first;
            elmnt != end();
            ++elmnt
        )
        {
            scalar de = mag(elmnt().position() - ip);
            scalar me = elmnt().m();
            mTot += me;

            n++;

            label i = 0;
            bool found = false;

            // insert the parcel in the correct place
            // and move the others 
            while ( ( i < n-1 ) && ( !found ) ) 
            {
                if (de < dist[i])
                {
                    found = true;
                    for(label j=n; j>i; j--)
                    {
                        m[j]    = m[j-1];
                        dist[j] = dist[j-1];
                    }
                    m[i]    = me;
                    dist[i] = de;
                }
                i++;
            }

            if (!found)
            {
                m[n]    = me;
                dist[n] = de;
            }
        }
    }

    reduce(mTot, sumOp<scalar>());

    if (Np > 1)
    {
        scalar mLimit = prc*mTot;
        scalar mOff = (1.0 - prc)*mTot;

        // 'prc' is large enough that the parcel most far
        // away will be used, no need to loop...
        if (mLimit > mTot - m[Np-1])
        {
            d = dist[Np-1];
        }
        else
        {
            scalar mOffSum = 0.0;
            label i = Np;

            while ((mOffSum < mOff) && (i>0))
            {
                i--;
                mOffSum += m[i];
            }
            d = dist[i];
        }

    }
    else
    {
        if (Np > 0)
        {
            spray::const_iterator elmnt = begin();
            d = mag(elmnt().position() - ip);
        }
    }

    reduce(d, maxOp<scalar>());

    return d;
}


scalar spray::smd() const
{
    scalar numerator = 0.0, denominator = VSMALL;

    for
    (
        spray::const_iterator elmnt = begin();
        elmnt != end();
        ++elmnt
    )
    {
        label celli = elmnt().cell();
        scalar Pc = p()[celli];
        scalar T = elmnt().T();
        scalar rho = fuels_->rho(Pc, T, elmnt().X());

        scalar tmp = elmnt().N(rho)*pow(elmnt().d(), 2.0);
        numerator += tmp*elmnt().d();
        denominator += tmp;
    }

    reduce(numerator, sumOp<scalar>());
    reduce(denominator, sumOp<scalar>());

    return numerator/denominator;
}


scalar spray::maxD() const
{
    scalar maxD = 0.0;

    for
    (
        spray::const_iterator elmnt = begin();
        elmnt != end();
        ++elmnt
    )
    {
        maxD = max(maxD, elmnt().d());
    }

    reduce(maxD, maxOp<scalar>());

    return maxD;
}


void spray::calculateAmbientPressure()
{
    ambientPressure_ = p_.average().value();
}


void spray::calculateAmbientTemperature()
{
    ambientTemperature_ = T_.average().value();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
