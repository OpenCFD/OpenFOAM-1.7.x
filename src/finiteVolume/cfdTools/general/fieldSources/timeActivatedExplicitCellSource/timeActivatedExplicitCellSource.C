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

#include "timeActivatedExplicitCellSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::timeActivatedExplicitCellSource::updateCellSet()
{
    cellSelector_->applyToSet(topoSetSource::NEW, selectedCellSet_);

    Info<< "    " << name_ << ": selected "
        << returnReduce(selectedCellSet_.size(), sumOp<label>())
        << " cells" << nl << endl;

    V_ = scalarField(selectedCellSet_.size(), 1.0);
    if (volumeType_ == vtAbsolute)
    {
        label i = 0;
        forAllConstIter(cellSet, selectedCellSet_, iter)
        {
            V_[i++] = mesh_.V()[iter.key()];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeActivatedExplicitCellSource::timeActivatedExplicitCellSource
(
    const word& name,
    const fvMesh& mesh,
    const dimensionSet& dims
)
:
    timeActivatedExplicitSource(name, mesh, dims),
    onValue_(readScalar(lookup("onValue"))),
    offValue_(readScalar(lookup("offValue"))),
    V_(0),
    cellSource_(lookup("cellSource")),
    cellSelector_
    (
        topoSetSource::New
        (
            cellSource_,
            mesh,
            subDict(cellSource_ + "Coeffs")
        )
    ),
    selectedCellSet_
    (
        mesh,
        name + "SourceCellSet",
        mesh.nCells()/10 + 1  // Reasonable size estimate.
    )
{
    // Create the cell set
    updateCellSet();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::timeActivatedExplicitCellSource::onValue() const
{
    return onValue_;
}


Foam::scalar Foam::timeActivatedExplicitCellSource::offValue() const
{
    return offValue_;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::timeActivatedExplicitCellSource::Su()
{
    tmp<DimensionedField<scalar, volMesh> > tSource
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                name_ + "Su",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensions_, 0.0)
        )
    );

    scalar value = offValue_;
    if
    (
        active_
     && (runTime_.time().value() >= timeStart_)
     && (runTime_.time().value() <= timeStart_ + duration_)
    )
    {
        // Update the cell set if the mesh is changing
        if (mesh_.changing())
        {
            updateCellSet();
        }

        value = onValue_;
    }

    DimensionedField<scalar, volMesh>& sourceField = tSource();

    label i = 0;
    forAllConstIter(cellSet, selectedCellSet_, iter)
    {
        sourceField[iter.key()] = value/V_[i++];
    }

    return tSource;
}


bool Foam::timeActivatedExplicitCellSource::read()
{
    if (timeActivatedExplicitSource::read())
    {
        lookup("onValue") >> onValue_;
        lookup("offValue") >> offValue_;
        lookup("cellSource") >> cellSource_;
        cellSelector_ =
            topoSetSource::New
            (
                cellSource_,
                mesh_,
                subDict(cellSource_ + "Coeffs")
            );
        updateCellSet();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
