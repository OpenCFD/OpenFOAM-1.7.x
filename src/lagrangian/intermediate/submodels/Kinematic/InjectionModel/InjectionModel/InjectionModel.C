/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "InjectionModel.H"
#include "mathematicalConstants.H"
#include "meshTools.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectionModel<CloudType>::readProps()
{
    IOobject propsDictHeader
    (
        "injectionProperties",
        owner_.db().time().timeName(),
        "uniform"/cloud::prefix/owner_.name(),
        owner_.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (propsDictHeader.headerOk())
    {
        const IOdictionary propsDict(propsDictHeader);

        propsDict.readIfPresent("massInjected", massInjected_);
        propsDict.readIfPresent("nInjections", nInjections_);
        propsDict.readIfPresent("parcelsAddedTotal", parcelsAddedTotal_);
        propsDict.readIfPresent("timeStep0", timeStep0_);
    }
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::writeProps()
{
    if (owner_.db().time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                "injectionProperties",
                owner_.db().time().timeName(),
                "uniform"/cloud::prefix/owner_.name(),
                owner_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        propsDict.add("massInjected", massInjected_);
        propsDict.add("nInjections", nInjections_);
        propsDict.add("parcelsAddedTotal", parcelsAddedTotal_);
        propsDict.add("timeStep0", timeStep0_);

        propsDict.regIOobject::write();
    }
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::prepareForNextTimeStep
(
    const scalar time,
    label& newParcels,
    scalar& newVolume
)
{
    // Initialise values
    newParcels = 0;
    newVolume = 0.0;

    // Return if not started injection event
    if (time < SOI_)
    {
        timeStep0_ = time;
        return;
    }

    // Make times relative to SOI
    scalar t0 = timeStep0_ - SOI_;
    scalar t1 = time - SOI_;

    // Number of parcels to inject
    newParcels = parcelsToInject(t0, t1);

    // Volume of parcels to inject
    newVolume = volumeToInject(t0, t1);

    // Hold previous time if no parcels, but non-zero volume fraction
    if ((newParcels == 0) && (newVolume > 0.0))
    {
        // hold value of timeStep0_
    }
    else
    {
        // advance value of timeStep0_
        timeStep0_ = time;
    }
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::findCellAtPosition
(
    label& cellI,
    vector& position
)
{
    const vector p0 = position;

    bool foundCell = false;

    cellI = owner_.mesh().findCell(position);

    if (cellI >= 0)
    {
        const vector& C = owner_.mesh().C()[cellI];
        position += SMALL*(C - position);

        foundCell = owner_.mesh().pointInCell(position, cellI);
    }
    reduce(foundCell, orOp<bool>());

    // Last chance - find nearest cell and try that one
    // - the point is probably on an edge
    if (!foundCell)
    {
        cellI = owner_.mesh().findNearestCell(position);

        if (cellI >= 0)
        {
            const vector& C = owner_.mesh().C()[cellI];
            position += SMALL*(C - position);

            foundCell = owner_.mesh().pointInCell(position, cellI);
        }
        reduce(foundCell, orOp<bool>());
    }

    if (!foundCell)
    {
        FatalErrorIn
        (
            "Foam::InjectionModel<CloudType>::findCellAtPosition"
            "("
                "label&, "
                "vector&"
            ")"
        )<< "Cannot find parcel injection cell. "
         << "Parcel position = " << p0 << nl
         << abort(FatalError);
    }
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::setNumberOfParticles
(
    const label parcels,
    const scalar volume,
    const scalar diameter,
    const scalar rho
)
{
    scalar nP = 0.0;
    switch (parcelBasis_)
    {
        case pbMass:
        {
            nP = volume/volumeTotal_
                *massTotal_/rho
               /(parcels*mathematicalConstant::pi/6.0*pow3(diameter));
            break;
        }
        case pbNumber:
        {
            nP = massTotal_/(rho*volumeTotal_*parcels);
            break;
        }
        default:
        {
            nP = 0.0;
            FatalErrorIn
            (
                "Foam::scalar "
                "Foam::InjectionModel<CloudType>::setNumberOfParticles"
                "("
                "    const label, "
                "    const scalar, "
                "    const scalar, "
                "    const scalar, "
                "    const scalar"
                ")"
            )<< "Unknown parcelBasis type" << nl
             << exit(FatalError);
        }
    }

    return nP;
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::postInjectCheck(const label parcelsAdded)
{
    if (parcelsAdded > 0)
    {
        Pout<< nl
            << "--> Cloud: " << owner_.name() << nl
            << "    Added " << parcelsAdded
            << " new parcels" << nl << endl;
    }

    // Increment total number of parcels added
    parcelsAddedTotal_ += parcelsAdded;

    // Update time for start of next injection
    time0_ = owner_.db().time().value();

    // Increment number of injections
    nInjections_++;

    // Write current state to properties file
    writeProps();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    SOI_(0.0),
    volumeTotal_(0.0),
    massTotal_(0.0),
    massInjected_(0.0),
    nInjections_(0),
    parcelsAddedTotal_(0),
    parcelBasis_(pbNumber),
    time0_(0.0),
    timeStep0_(0.0)
{
    readProps();
}


template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    SOI_(readScalar(coeffDict_.lookup("SOI"))),
    volumeTotal_(0.0),
    massTotal_(dimensionedScalar(coeffDict_.lookup("massTotal")).value()),
    massInjected_(0.0),
    nInjections_(0),
    parcelsAddedTotal_(0),
    parcelBasis_(pbNumber),
    time0_(owner.db().time().value()),
    timeStep0_(0.0)
{
    // Provide some info
    // - also serves to initialise mesh dimensions - needed for parallel runs
    //   due to lazy evaluation of valid mesh dimensions
    Info<< "    Constructing " << owner.mesh().nGeometricD() << "-D injection"
        << endl;

    word parcelBasisType = coeffDict_.lookup("parcelBasisType");
    if (parcelBasisType == "mass")
    {
        parcelBasis_ = pbMass;
    }
    else if (parcelBasisType == "number")
    {
        parcelBasis_ = pbNumber;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::InjectionModel<CloudType>::InjectionModel"
            "("
                "const dictionary&, "
                "CloudType&, "
                "const word&"
            ")"
        )<< "parcelBasisType must be either 'number' or 'mass'" << nl
         << exit(FatalError);
    }

    readProps();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::~InjectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class TrackData>
void Foam::InjectionModel<CloudType>::inject(TrackData& td)
{
    if (!active())
    {
        return;
    }

    const scalar time = owner_.db().time().value();
    const scalar carrierDt = owner_.db().time().deltaT().value();
    const polyMesh& mesh = owner_.mesh();

    // Prepare for next time step
    label newParcels = 0;
    scalar newVolume = 0.0;
    prepareForNextTimeStep(time, newParcels, newVolume);

    // Return if no parcels are required
    if (newParcels == 0)
    {
        postInjectCheck(0);
        return;
    }

    // Duration of injection period during this timestep
    const scalar deltaT =
        max(0.0, min(carrierDt, min(time - SOI_, timeEnd() - time0_)));

    // Pad injection time if injection starts during this timestep
    const scalar padTime = max(0.0, SOI_ - time0_);

    // Introduce new parcels linearly across carrier phase timestep
    label parcelsAdded = 0;
    for (label parcelI=0; parcelI<newParcels; parcelI++)
    {
        if (validInjection(parcelI))
        {
            // Calculate the pseudo time of injection for parcel 'parcelI'
            scalar timeInj = time0_ + padTime + deltaT*parcelI/newParcels;

            // Determine the injection position and owner cell
            label cellI = -1;
            vector pos = vector::zero;
            setPositionAndCell(parcelI, newParcels, timeInj, pos, cellI);

            if (cellI > -1)
            {
                // Lagrangian timestep
                scalar dt = time - timeInj;

                // Apply corrections to position for 2-D cases
                meshTools::constrainToMeshCentre(mesh, pos);

                // Create a new parcel
                parcelType* pPtr = new parcelType(td.cloud(), pos, cellI);

                // Assign new parcel properties in injection model
                setProperties(parcelI, newParcels, timeInj, *pPtr);

                // Check new parcel properties
                td.cloud().checkParcelProperties(*pPtr, dt, fullyDescribed());

                // Apply correction to velocity for 2-D cases
                meshTools::constrainDirection
                (
                    mesh,
                    mesh.solutionD(),
                    pPtr->U()
                );

                // Number of particles per parcel
                pPtr->nParticle() =
                    setNumberOfParticles
                    (
                        newParcels,
                        newVolume,
                        pPtr->d(),
                        pPtr->rho()
                    );

                // Add the new parcel
                td.cloud().addParticle(pPtr);

                massInjected_ += pPtr->nParticle()*pPtr->mass();
                parcelsAdded++;
            }
        }
    }

    postInjectCheck(parcelsAdded);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewInjectionModel.C"

// ************************************************************************* //
