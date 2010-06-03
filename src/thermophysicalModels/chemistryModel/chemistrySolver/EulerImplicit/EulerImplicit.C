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

#include "EulerImplicit.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::EulerImplicit<CompType, ThermoType>::EulerImplicit
(
    ODEChemistryModel<CompType, ThermoType>& model,
    const word& modelName
)
:
    chemistrySolver<CompType, ThermoType>(model, modelName),
    coeffsDict_(model.subDict(modelName + "Coeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    equil_(coeffsDict_.lookup("equilibriumRateLimiter"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::EulerImplicit<CompType, ThermoType>::~EulerImplicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::scalar Foam::EulerImplicit<CompType, ThermoType>::solve
(
    scalarField &c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    label nSpecie = this->model_.nSpecie();
    simpleMatrix<scalar> RR(nSpecie);

    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    for (label i=0; i<nSpecie; i++)
    {
        RR.source()[i] = c[i]/dt;
    }

    for (label i=0; i<this->model_.reactions().size(); i++)
    {
        const Reaction<ThermoType>& R = this->model_.reactions()[i];

        scalar omegai = this->model_.omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        scalar corr = 1.0;
        if (equil_)
        {
            if (omegai<0.0)
            {
                corr = 1.0/(1.0 + pr*dt);
            }
            else
            {
                corr = 1.0/(1.0 + pf*dt);
            }
        }

        for (label s=0; s<R.lhs().size(); s++)
        {
            label si = R.lhs()[s].index;
            scalar sl = R.lhs()[s].stoichCoeff;
            RR[si][rRef] -= sl*pr*corr;
            RR[si][lRef] += sl*pf*corr;
        }

        for (label s=0; s<R.rhs().size(); s++)
        {
            label si = R.rhs()[s].index;
            scalar sr = R.rhs()[s].stoichCoeff;
            RR[si][lRef] -= sr*pf*corr;
            RR[si][rRef] += sr*pr*corr;
        }

    } // end for(label i...


    for (label i=0; i<nSpecie; i++)
    {
        RR[i][i] += 1.0/dt;
    }

    c = RR.LUsolve();
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    // estimate the next time step
    scalar tMin = GREAT;
    label nEqns = this->model_.nEqns();
    scalarField c1(nEqns, 0.0);

    for (label i=0; i<nSpecie; i++)
    {
        c1[i] = c[i];
    }
    c1[nSpecie] = T;
    c1[nSpecie+1] = p;

    scalarField dcdt(nEqns, 0.0);
    this->model_.derivatives(0.0, c1, dcdt);

    scalar sumC = sum(c);

    for (label i=0; i<nSpecie; i++)
    {
        scalar d = dcdt[i];
        if (d < -SMALL)
        {
            tMin = min(tMin, -(c[i] + SMALL)/d);
        }
        else
        {
            d = max(d, SMALL);
            scalar cm = max(sumC - c[i], 1.0e-5);
            tMin = min(tMin, cm/d);
        }
    }

    return cTauChem_*tMin;
}


// ************************************************************************* //
