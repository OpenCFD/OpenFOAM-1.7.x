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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::word& Foam::DsmcCloud<ParcelType>::cloudName() const
{
    return cloudName_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::DsmcCloud<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::DsmcCloud<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const Foam::List<Foam::word>&
Foam::DsmcCloud<ParcelType>::typeIdList() const
{
    return typeIdList_;
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline const Foam::List<Foam::DynamicList<ParcelType*> >&
Foam::DsmcCloud<ParcelType>::cellOccupancy() const
{
    return cellOccupancy_;
}


template<class ParcelType>
inline Foam::volScalarField& Foam::DsmcCloud<ParcelType>::sigmaTcRMax()
{
    return sigmaTcRMax_;
}


template<class ParcelType>
inline Foam::scalarField&
Foam::DsmcCloud<ParcelType>::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


template<class ParcelType>
inline const Foam::List<typename ParcelType::constantProperties>&
Foam::DsmcCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::DsmcCloud<ParcelType>::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorIn("Foam::DsmcCloud<ParcelType>::constProps(label typeId)")
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}


template<class ParcelType>
inline Foam::Random& Foam::DsmcCloud<ParcelType>::rndGen()
{
    return rndGen_;
}


template<class ParcelType>
inline void Foam::DsmcCloud<ParcelType>::storeDeltaT()
{
    cachedDeltaT_ = mesh().time().deltaT().value();
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::cachedDeltaT() const
{
    return cachedDeltaT_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::DsmcCloud<ParcelType>::q() const
{
    return q_;
}


template<class ParcelType>
inline Foam::volScalarField& Foam::DsmcCloud<ParcelType>::q()
{
    return q_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::DsmcCloud<ParcelType>::fD() const
{
    return fD_;
}


template<class ParcelType>
inline Foam::volVectorField& Foam::DsmcCloud<ParcelType>::fD()
{
    return fD_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::DsmcCloud<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::DsmcCloud<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline const Foam::BinaryCollisionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::binaryCollision() const
{
    return binaryCollisionModel_;
}


template<class ParcelType>
inline Foam::BinaryCollisionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::binaryCollision()
{
    return binaryCollisionModel_();
}


template<class ParcelType>
inline const Foam::WallInteractionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::wallInteraction() const
{
    return wallInteractionModel_;
}


template<class ParcelType>
inline Foam::WallInteractionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::wallInteraction()
{
    return wallInteractionModel_();
}


template<class ParcelType>
inline const Foam::InflowBoundaryModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::inflowBoundary() const
{
    return inflowBoundaryModel_;
}


template<class ParcelType>
inline Foam::InflowBoundaryModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::inflowBoundary()
{
    return inflowBoundaryModel_();
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        sysMass += cP.mass();
    }

    return nParticle_*sysMass;
}


template<class ParcelType>
inline Foam::vector Foam::DsmcCloud<ParcelType>::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearMomentum += cP.mass()*p.U();
    }

    return nParticle_*linearMomentum;
}


template<class ParcelType>
inline Foam::scalar
Foam::DsmcCloud<ParcelType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearKineticEnergy += 0.5*cP.mass()*(p.U() & p.U());
    }

    return nParticle_*linearKineticEnergy;
}


template<class ParcelType>
inline Foam::scalar
Foam::DsmcCloud<ParcelType>::internalEnergyOfSystem() const
{
    scalar internalEnergy = 0.0;

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        internalEnergy += p.Ei();
    }

    return nParticle_*internalEnergy;
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::maxwellianAverageSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return 2.0*sqrt(2.0*kb*temperature/(mathematicalConstant::pi*mass));
}


template<class ParcelType>
inline Foam::scalarField Foam::DsmcCloud<ParcelType>::maxwellianAverageSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return 2.0*sqrt(2.0*kb*temperature/(mathematicalConstant::pi*mass));
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::maxwellianRMSSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(3.0*kb*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField Foam::DsmcCloud<ParcelType>::maxwellianRMSSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return sqrt(3.0*kb*temperature/mass);
}


template<class ParcelType>
inline Foam::scalar
Foam::DsmcCloud<ParcelType>::maxwellianMostProbableSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(2.0*kb*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField
Foam::DsmcCloud<ParcelType>::maxwellianMostProbableSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return sqrt(2.0*kb*temperature/mass);
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::rhoN() const
{
    tmp<volScalarField> trhoN
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "rhoN",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, -3, 0, 0, 0), VSMALL)
        )
    );

    scalarField& rhoN = trhoN().internalField();
    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        rhoN[cellI]++;
    }

    rhoN *= nParticle_/mesh().cellVolumes();

    return trhoN;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::rhoM() const
{
    tmp<volScalarField> trhoM
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "rhoM",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0), VSMALL)
        )
    );

    scalarField& rhoM = trhoM().internalField();
    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        rhoM[cellI] += constProps(p.typeId()).mass();
    }

    rhoM *= nParticle_/mesh().cellVolumes();

    return trhoM;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::dsmcRhoN() const
{
    tmp<volScalarField> tdsmcRhoN
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "dsmcRhoN",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, -3, 0, 0, 0), 0.0)
        )
    );

    scalarField& dsmcRhoN = tdsmcRhoN().internalField();
    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        dsmcRhoN[cellI]++;
    }

    return tdsmcRhoN;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volVectorField>
Foam::DsmcCloud<ParcelType>::momentum() const
{
    tmp<volVectorField> tmomentum
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + "momentum",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector
            (
                "zero",
                dimensionSet(1, -2, -1, 0, 0),
                vector::zero
            )
        )
    );

    vectorField& momentum = tmomentum().internalField();
    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        momentum[cellI] += constProps(p.typeId()).mass()*p.U();
    }

    momentum *= nParticle_/mesh().cellVolumes();

    return tmomentum;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::linearKE() const
{
    tmp<volScalarField> tlinearKE
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "linearKE",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
        )
    );

    scalarField& linearKE = tlinearKE().internalField();
    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        linearKE[cellI] += 0.5*constProps(p.typeId()).mass()*(p.U() & p.U());
    }

    linearKE *= nParticle_/mesh().cellVolumes();

    return tlinearKE;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::internalE() const
{
    tmp<volScalarField> tinternalE
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "internalE",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
        )
    );

    scalarField& internalE = tinternalE().internalField();
    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        internalE[cellI] += p.Ei();
    }

    internalE *= nParticle_/mesh().cellVolumes();

    return tinternalE;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::iDof() const
{
    tmp<volScalarField> tiDof
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "iDof",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, -3, 0, 0, 0), VSMALL)
        )
    );

    scalarField& iDof = tiDof().internalField();

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        iDof[cellI] += constProps(p.typeId()).internalDegreesOfFreedom();
    }

    iDof *= nParticle_/mesh().cellVolumes();

    return tiDof;
}


template<class ParcelType>
inline void Foam::DsmcCloud<ParcelType>::clear()
{
    return IDLList<ParcelType>::clear();
}


// ************************************************************************* //
