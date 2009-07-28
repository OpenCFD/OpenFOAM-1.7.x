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

#include "timeActivatedExplicitMulticomponentPointSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label
Foam::timeActivatedExplicitMulticomponentPointSource::carrierFieldId
(
    const word& fieldName
)
{
    forAll(carrierFields_, fieldI)
    {
        if (carrierFields_[fieldI].name() == fieldName)
        {
            return fieldI;
        }
    }

    return -1;
}


void Foam::timeActivatedExplicitMulticomponentPointSource::updateAddressing()
{
    forAll(pointSources_, sourceI)
    {
        const pointSourceProperties& psp = pointSources_[sourceI];
        bool foundCell = false;
        label cid = mesh_.findCell(psp.location());
        if (cid >= 0)
        {
            foundCell = mesh_.pointInCell(psp.location(), cid);
        }
        reduce(foundCell, orOp<bool>());
        if (!foundCell)
        {
            label cid = mesh_.findNearestCell(psp.location());
            if (cid >= 0)
            {
                foundCell = mesh_.pointInCell(psp.location(), cid);
            }
        }
        reduce(foundCell, orOp<bool>());

        if (!foundCell)
        {
            FatalErrorIn
            (
                "timeActivatedExplicitMulticomponentPointSource::"
                "updateAddressing()"
            )   << "Unable to find location " << psp.location() << " in mesh "
                << "for source " << psp.name() << nl
                << exit(FatalError);
        }
        else
        {
            cellOwners_[sourceI] = cid;
        }

        fieldIds_[sourceI].setSize(psp.fieldData().size());
        forAll(psp.fieldData(), fieldI)
        {
            const word& fieldName = psp.fieldData()[fieldI].first();
            label cfid = carrierFieldId(fieldName);
            if (cfid < 0)
            {
                FatalErrorIn
                (
                    "timeActivatedExplicitMulticomponentPointSource::"
                    "updateAddressing()"
                )   << "Unable to find field " << fieldName << " in carrier "
                    << "fields for source " << psp.name() << nl
                    << exit(FatalError);
            }
            else
            {
                fieldIds_[sourceI][fieldI] = cfid;
            }
       }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeActivatedExplicitMulticomponentPointSource::
timeActivatedExplicitMulticomponentPointSource
(
    const word& name,
    const fvMesh& mesh,
    const PtrList<volScalarField>& carrierFields,
    const dimensionSet& dims
)
:
    IOdictionary
    (
        IOobject
        (
            name + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    mesh_(mesh),
    runTime_(mesh.time()),
    dimensions_(dims),
    carrierFields_(carrierFields),
    active_(lookup("active")),
    pointSources_(lookup("pointSources")),
    cellOwners_(pointSources_.size()),
    fieldIds_(pointSources_.size())
{
    // Initialise the field addressing
    updateAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::timeActivatedExplicitMulticomponentPointSource::Su
(
    const label fieldI
)
{
    if (mesh_.changing())
    {
        updateAddressing();
    }

    tmp<DimensionedField<scalar, volMesh> > tSource
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name_ + carrierFields_[fieldI].name() + "Su",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensions_, 0.0)
        )
    );

    if (active_)
    {
        DimensionedField<scalar, volMesh>& sourceField = tSource();

        const scalarField& V = mesh_.V();
        const scalar dt = runTime_.deltaT().value();

        forAll(pointSources_, sourceI)
        {
            const pointSourceProperties& psp = pointSources_[sourceI];

            forAll(fieldIds_[sourceI], i)
            {
                if
                (
                    fieldIds_[sourceI][i] == fieldI
                && (runTime_.time().value() >= psp.timeStart())
                && (runTime_.time().value() <= psp.timeEnd())
                )
                {
                    const label cid = cellOwners_[sourceI];
                    if (cid >= 0)
                    {
                        sourceField[cid] +=
                            dt*psp.fieldData()[i].second()/V[cid];
                    }
                }
            }
        }
    }

    return tSource;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::timeActivatedExplicitMulticomponentPointSource::Su()
{
    if (mesh_.changing())
    {
        updateAddressing();
    }

    tmp<DimensionedField<scalar, volMesh> > tSource
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name_ + "TotalSu",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensions_, 0.0)
        )
    );

    if (active_)
    {
        DimensionedField<scalar, volMesh>& sourceField = tSource();

        const scalarField& V = mesh_.V();
        const scalar dt = runTime_.deltaT().value();

        forAll(pointSources_, sourceI)
        {
            const pointSourceProperties& psp = pointSources_[sourceI];

            forAll(fieldIds_[sourceI], i)
            {
                if
                (
                   (runTime_.time().value() >= psp.timeStart())
                && (runTime_.time().value() <= psp.timeEnd())
                )
                {
                    const label cid = cellOwners_[sourceI];
                    if (cid >= 0)
                    {
                        sourceField[cid] +=
                            dt*psp.fieldData()[i].second()/V[cid];
                    }
                }
            }
        }
    }

    return tSource;
}


bool Foam::timeActivatedExplicitMulticomponentPointSource::read()
{
    if (regIOobject::read())
    {
        lookup("active") >> active_;
        lookup("pointSources") >> pointSources_;

        cellOwners_.setSize(pointSources_.size());

        updateAddressing();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
