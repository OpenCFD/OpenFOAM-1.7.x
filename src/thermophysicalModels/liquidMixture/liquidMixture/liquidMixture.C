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

#include "liquidMixture.H"
#include "dictionary.H"
#include "specie.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::liquidMixture::TrMax = 0.999;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidMixture::liquidMixture
(
    const dictionary& thermophysicalProperties
)
:
    components_(thermophysicalProperties.lookup("liquidComponents")),
    properties_(components_.size())
{
    // use sub-dictionary "liquidProperties" if possible to avoid
    // collisions with identically named gas-phase entries
    // (eg, H2O liquid vs. gas)
    forAll(components_, i)
    {
        const dictionary* subDictPtr = thermophysicalProperties.subDictPtr
        (
            "liquidProperties"
        );

        if (subDictPtr)
        {
            properties_.set
            (
                i,
                liquid::New(subDictPtr->lookup(components_[i]))
            );
        }
        else
        {
            properties_.set
            (
                i,
                liquid::New(thermophysicalProperties.lookup(components_[i]))
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liquidMixture> Foam::liquidMixture::New
(
    const dictionary& thermophysicalProperties
)
{
    return autoPtr<liquidMixture>(new liquidMixture(thermophysicalProperties));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::liquidMixture::Tc
(
    const scalarField& x
) const
{

    scalar vTc = 0.0;
    scalar vc = 0.0;

    forAll(properties_, i)
    {
        scalar x1 = x[i]*properties_[i].Vc();
        vc += x1;
        vTc += x1*properties_[i].Tc();
    }

    return vTc/vc;
}


Foam::scalar Foam::liquidMixture::Tpc
(
    const scalarField& x
) const
{
    scalar Tpc = 0.0;
    forAll(properties_, i)
    {
        Tpc += x[i]*properties_[i].Tc();
    }

    return Tpc;
}


Foam::scalar Foam::liquidMixture::Ppc
(
    const scalarField& x
) const
{
    scalar Vc = 0.0;
    scalar Zc = 0.0;
    forAll(properties_, i)
    {
        Vc += x[i]*properties_[i].Vc();
        Zc += x[i]*properties_[i].Zc();
    }

    return specie::RR*Zc*Tpc(x)/Vc;
}


Foam::scalar Foam::liquidMixture::omega
(
    const scalarField& x
) const
{
    scalar omega = 0.0;
    forAll(properties_, i)
    {
        omega += x[i]*properties_[i].omega();
    }

    return omega;
}


Foam::scalarField Foam::liquidMixture::Xs
(
    const scalar p,
    const scalar Tg,
    const scalar Tl,
    const scalarField& xg,
    const scalarField& xl
) const
{
    scalarField xs(xl.size(), 0.0);

    // Raoult's Law
    forAll(xs, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), Tl);
        xs[i] = properties_[i].pv(p, Ti)*xl[i]/p;
    }
    return xs;
}


Foam::scalar Foam::liquidMixture::W
(
    const scalarField& x
) const
{
    scalar W = 0.0;
    forAll(properties_, i)
    {
        W += x[i]*properties_[i].W();
    }

    return W;
}


Foam::scalarField Foam::liquidMixture::Y
(
    const scalarField& X
) const
{
    scalarField Y(X/W(X));

    forAll(Y, i)
    {
        Y[i] *= properties_[i].W();
    }

    return Y;
}


Foam::scalarField Foam::liquidMixture::X
(
    const scalarField& Y
) const
{
    scalarField X(Y.size());
    scalar Winv = 0.0;
    forAll(X, i)
    {
        Winv += Y[i]/properties_[i].W();
        X[i] = Y[i]/properties_[i].W();
    }

    return X/Winv;
}


Foam::scalar Foam::liquidMixture::rho
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar v = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            scalar rho = SMALL + properties_[i].rho(p, Ti);
            v += x[i]*properties_[i].W()/rho;
        }
    }

    return W(x)/v;
}


Foam::scalar Foam::liquidMixture::pv
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar pv = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            pv += x[i]*properties_[i].pv(p, Ti)*properties_[i].W();
        }
    }

    return pv/W(x);
}


Foam::scalar Foam::liquidMixture::hl
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar hl = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            hl += x[i]*properties_[i].hl(p, Ti)*properties_[i].W();
        }
    }

    return hl/W(x);
}


Foam::scalar Foam::liquidMixture::cp
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar cp = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            cp += x[i]*properties_[i].cp(p, Ti)*properties_[i].W();
        }
    }

    return cp/W(x);
}


Foam::scalar Foam::liquidMixture::sigma
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    // sigma is based on surface mole fractions
    // which is estimated from Raoult's Law
    scalar sigma = 0.0;
    scalarField Xs(x.size(), 0.0);
    scalar XsSum = 0.0;
    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);
        scalar Pvs = properties_[i].pv(p, Ti);
        scalar xs = x[i]*Pvs/p;
        XsSum += xs;
        Xs[i] = xs;
    }

    forAll(properties_, i)
    {
        if (Xs[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            sigma += (Xs[i]/XsSum)*properties_[i].sigma(p, Ti);
        }
    }

    return sigma;
}


Foam::scalar Foam::liquidMixture::mu
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar mu = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            mu += x[i]*log(properties_[i].mu(p, Ti));
        }
    }

    return exp(mu);
}


Foam::scalar Foam::liquidMixture::K
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    // calculate superficial volume fractions phii
    scalarField phii(x.size(), 0.0);
    scalar pSum = 0.0;

    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);

        scalar Vi = properties_[i].W()/properties_[i].rho(p, Ti);
        phii[i] = x[i]*Vi;
        pSum += phii[i];
    }

    forAll(phii, i)
    {
        phii[i] /= pSum;
    }

    scalar K = 0.0;

    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);

        forAll(properties_, j)
        {
            scalar Tj = min(TrMax*properties_[j].Tc(), T);

            scalar Kij =
                2.0
               /(
                    1.0/properties_[i].K(p, Ti)
                  + 1.0/properties_[j].K(p, Tj)
                );
            K += phii[i]*phii[j]*Kij;
        }
    }

    return K;
}


Foam::scalar Foam::liquidMixture::D
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    // Blanc's law
    scalar Dinv = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            Dinv += x[i]/properties_[i].D(p, Ti);
        }
    }

    return 1.0/Dinv;
}


// ************************************************************************* //
