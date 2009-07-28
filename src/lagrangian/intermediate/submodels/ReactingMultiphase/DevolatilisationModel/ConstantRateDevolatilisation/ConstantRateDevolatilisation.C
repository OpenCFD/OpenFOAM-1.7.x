/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "ConstantRateDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ConstantRateDevolatilisation<CloudType>::ConstantRateDevolatilisation
(
    const dictionary& dict,
    CloudType& owner
)
:
    DevolatilisationModel<CloudType>(dict, owner, typeName),
    A0_(dimensionedScalar(this->coeffDict().lookup("A0")).value()),
    volatileResidualCoeff_
    (
        readScalar(this->coeffDict().lookup("volatileResidualCoeff"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ConstantRateDevolatilisation<CloudType>::~ConstantRateDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ConstantRateDevolatilisation<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::ConstantRateDevolatilisation<CloudType>::calculate
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
    const scalar massVolatile0 = YVolatile0*mass0;
    const scalar massVolatile  = YVolatile*mass;

    if (massVolatile <= volatileResidualCoeff_*massVolatile0)
    {
        canCombust = true;
    }

    // Volatile devolatilisation from particle to carrier gas phase
    const scalar dMass = min(dt*A0_*massVolatile0, massVolatile);

    return dMass;
}


// ************************************************************************* //
