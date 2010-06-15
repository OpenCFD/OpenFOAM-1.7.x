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

#include "ode.H"
#include "ODEChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ode<CompType, ThermoType>::ode
(
    ODEChemistryModel<CompType, ThermoType>& model,
    const word& modelName
)
:
    chemistrySolver<CompType, ThermoType>(model, modelName),
    coeffsDict_(model.subDict(modelName + "Coeffs")),
    solverName_(coeffsDict_.lookup("ODESolver")),
    odeSolver_(ODESolver::New(solverName_, model)),
    eps_(readScalar(coeffsDict_.lookup("eps"))),
    scale_(readScalar(coeffsDict_.lookup("scale")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ode<CompType, ThermoType>::~ode()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::scalar Foam::ode<CompType, ThermoType>::solve
(
    scalarField& c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    label nSpecie = this->model_.nSpecie();
    scalarField c1(this->model_.nEqns(), 0.0);

    // copy the concentration, T and P to the total solve-vector
    for (label i=0; i<nSpecie; i++)
    {
        c1[i] = c[i];
    }
    c1[nSpecie] = T;
    c1[nSpecie+1] = p;

    scalar dtEst = dt;

    odeSolver_->solve
    (
        this->model_,
        t0,
        t0 + dt,
        c1,
        eps_,
        dtEst
    );

    for (label i=0; i<c.size(); i++)
    {
        c[i] = max(0.0, c1[i]);
    }

    return dtEst;
}


// ************************************************************************* //
