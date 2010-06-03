/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "SingleKineticRateDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SingleKineticRateDevolatilisation<CloudType>::
SingleKineticRateDevolatilisation
(
    const dictionary& dict,
    CloudType& owner
)
:
    DevolatilisationModel<CloudType>(dict, owner, typeName),
    A1_(dimensionedScalar(this->coeffDict().lookup("A1")).value()),
    E_(dimensionedScalar(this->coeffDict().lookup("E")).value()),
    volatileResidualCoeff_
    (
        readScalar(this->coeffDict().lookup("volatileResidualCoeff"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SingleKineticRateDevolatilisation<CloudType>::
~SingleKineticRateDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::SingleKineticRateDevolatilisation<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::SingleKineticRateDevolatilisation<CloudType>::calculate
(
    const scalar dt,
    const scalar mass0,
    const scalar mass,
    const scalar T,
    const scalar YVolatile0,
    const scalar YVolatile,
    bool& canCombust
) const
{
    const scalar massVolatile0 = YVolatile0*mass;
    const scalar massVolatile = YVolatile*mass;

    if (massVolatile <= volatileResidualCoeff_*massVolatile0)
    {
        canCombust = true;
    }

    // Kinetic rate
    const scalar kappa = A1_*exp(-E_/(specie::RR*T));

    // Volatile devolatilisation from particle to carrier gas phase
    const scalar dMass = min(dt*kappa*massVolatile, massVolatile);

    return dMass;
}


// ************************************************************************* //
