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

#include "pressureGradientExplicitSource.H"
#include "volFields.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pressureGradientExplicitSource::writeGradP() const
{
    // Only write on output time
    if (mesh_.time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                sourceName_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP_);
        propsDict.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureGradientExplicitSource::pressureGradientExplicitSource
(
    const word& sourceName,
    volVectorField& U
)
:
    sourceName_(sourceName),
    mesh_(U.mesh()),
    U_(U),
    dict_
    (
        IOobject
        (
            sourceName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    Ubar_(dict_.lookup("Ubar")),
    gradPini_(dict_.lookup("gradPini")),
    gradP_(gradPini_),
    flowDir_(Ubar_/mag(Ubar_)),
    cellSource_(dict_.lookup("cellSource")),
    cellSelector_
    (
        topoSetSource::New
        (
            cellSource_,
            mesh_,
            dict_.subDict(cellSource_ + "Coeffs")
        )
    ),
    selectedCellSet_
    (
        mesh_,
        sourceName_ + "CellSet",
        mesh_.nCells()/10 + 1  // Reasonable size estimate.
    )
{
    // Create the cell set
    cellSelector_->applyToSet
    (
        topoSetSource::NEW,
        selectedCellSet_
    );

    // Give some feedback
    Info<< "    Selected "
        << returnReduce(selectedCellSet_.size(), sumOp<label>())
        << " cells" << endl;

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timeName()/"uniform"/(sourceName_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP_;
    }

    Info<< "    Initial pressure gradient = " << gradP_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh> >
Foam::pressureGradientExplicitSource::Su() const
{
    tmp<DimensionedField<vector, volMesh> > tSource
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                sourceName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", gradP_.dimensions(), vector::zero)
        )
    );

    DimensionedField<vector, volMesh>& sourceField = tSource();

    forAllConstIter(cellSet, selectedCellSet_, iter)
    {
        label cellI = iter.key();

        sourceField[cellI] = flowDir_*gradP_.value();
    }

    return tSource;
}


void Foam::pressureGradientExplicitSource::update()
{
    const volScalarField& rUA =
        mesh_.lookupObject<volScalarField>("(1|A(" + U_.name() + "))");

    // Integrate flow variables over cell set
    scalar volTot = 0.0;
    scalar magUbarAve = 0.0;
    scalar rUAave = 0.0;
    forAllConstIter(cellSet, selectedCellSet_, iter)
    {
        label cellI = iter.key();

        scalar volCell = mesh_.V()[cellI];
        volTot += volCell;

        magUbarAve += (flowDir_ & U_[cellI])*volCell;
        rUAave += rUA[cellI]*volCell;
    }

    // Collect across all processors
    reduce(volTot, sumOp<scalar>());
    reduce(magUbarAve, sumOp<scalar>());
    reduce(rUAave, sumOp<scalar>());

    // Volume averages
    magUbarAve /= volTot;
    rUAave /= volTot;

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    scalar gradPplus = (mag(Ubar_) - magUbarAve)/rUAave;

    // Apply correction to velocity field
    forAllConstIter(cellSet, selectedCellSet_, iter)
    {
        label cellI = iter.key();
        U_[cellI] += flowDir_*rUA[cellI]*gradPplus;
    }

    // Update pressure gradient
    gradP_.value() += gradPplus;

    Info<< "Uncorrected Ubar = " << magUbarAve << tab
        << "Pressure gradient = " << gradP_.value() << endl;

    writeGradP();
}


// ************************************************************************* //
