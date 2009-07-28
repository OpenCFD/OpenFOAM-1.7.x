/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "LarsenBorgnakkeVariableHardSphere.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
Foam::scalar Foam::LarsenBorgnakkeVariableHardSphere<CloudType>::energyRatio
(
    scalar ChiA,
    scalar ChiB
)
{
    CloudType& cloud(this->owner());

    Random& rndGen(cloud.rndGen());

    scalar ChiAMinusOne = ChiA - 1;

    scalar ChiBMinusOne = ChiB - 1;

    if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
    {
        return rndGen.scalar01();
    }

    scalar energyRatio;

    scalar P;

    do
    {
        P = 0;

        energyRatio = rndGen.scalar01();

        if (ChiAMinusOne < SMALL)
        {
            P = 1.0 - pow(energyRatio, ChiB);
        }
        else if (ChiBMinusOne < SMALL)
        {
            P = 1.0 - pow(energyRatio, ChiA);
        }
        else
        {
            P =
                pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne,
                    ChiAMinusOne
                )
               *pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*(1 - energyRatio)
                    /ChiBMinusOne,
                    ChiBMinusOne
                );
        }
    } while (P < rndGen.scalar01());

    return energyRatio;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LarsenBorgnakkeVariableHardSphere<CloudType>::LarsenBorgnakkeVariableHardSphere
(
    const dictionary& dict,
    CloudType& cloud
)
:
    BinaryCollisionModel<CloudType>(dict, cloud, typeName),
    Tref_(readScalar(this->coeffDict().lookup("Tref"))),
    relaxationCollisionNumber_
    (
        readScalar(this->coeffDict().lookup("relaxationCollisionNumber"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LarsenBorgnakkeVariableHardSphere<CloudType>::
~LarsenBorgnakkeVariableHardSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template <class CloudType>
Foam::scalar Foam::LarsenBorgnakkeVariableHardSphere<CloudType>::sigmaTcR
(
    label typeIdP,
    label typeIdQ,
    const vector& UP,
    const vector& UQ
) const
{
    const CloudType& cloud(this->owner());

    scalar dPQ =
        0.5
       *(
            cloud.constProps(typeIdP).d()
          + cloud.constProps(typeIdQ).d()
        );

    scalar omegaPQ =
        0.5
       *(
            cloud.constProps(typeIdP).omega()
          + cloud.constProps(typeIdQ).omega()
        );

    scalar cR = mag(UP - UQ);

    if (cR < VSMALL)
    {
        return 0;
    }

    scalar mP = cloud.constProps(typeIdP).mass();

    scalar mQ = cloud.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    // calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    scalar sigmaTPQ =
        mathematicalConstant::pi*dPQ*dPQ
       *pow(2.0*CloudType::kb*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}


template <class CloudType>
void Foam::LarsenBorgnakkeVariableHardSphere<CloudType>::collide
(
    label typeIdP,
    label typeIdQ,
    vector& UP,
    vector& UQ,
    scalar& EiP,
    scalar& EiQ
)
{
    CloudType& cloud(this->owner());

    Random& rndGen(cloud.rndGen());

    scalar inverseCollisionNumber = 1/relaxationCollisionNumber_;

    // Larsen Borgnakke internal energy redistribution part.  Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // DSMC0R.FOR

    scalar preCollisionEiP = EiP;

    scalar preCollisionEiQ = EiQ;

    scalar iDofP = cloud.constProps(typeIdP).internalDegreesOfFreedom();

    scalar iDofQ = cloud.constProps(typeIdQ).internalDegreesOfFreedom();

    scalar omegaPQ =
        0.5
       *(
            cloud.constProps(typeIdP).omega()
          + cloud.constProps(typeIdQ).omega()
        );

    scalar mP = cloud.constProps(typeIdP).mass();

    scalar mQ = cloud.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

    scalar cRsqr = magSqr(UP - UQ);

    scalar availableEnergy = 0.5*mR*cRsqr;

    scalar ChiB = 2.5 - omegaPQ;

    if (iDofP > 0)
    {
        if (inverseCollisionNumber > rndGen.scalar01())
        {
            availableEnergy += preCollisionEiP;

            scalar ChiA = 0.5*iDofP;

            EiP = energyRatio(ChiA, ChiB)*availableEnergy;

            availableEnergy -= EiP;
        }
    }

    if (iDofQ > 0)
    {
        if (inverseCollisionNumber > rndGen.scalar01())
        {
            availableEnergy += preCollisionEiQ;

            // Change to general LB ratio calculation
            scalar energyRatio = 1.0 - pow(rndGen.scalar01(),(1.0/ChiB));

            EiQ = energyRatio*availableEnergy;

            availableEnergy -= EiQ;
        }
    }

    // Rescale the translational energy
    scalar cR = sqrt(2.0*availableEnergy/mR);

    // Variable Hard Sphere collision part

    scalar cosTheta = 2.0*rndGen.scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = 2.0*mathematicalConstant::pi*rndGen.scalar01();

    vector postCollisionRelU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
}


// ************************************************************************* //
