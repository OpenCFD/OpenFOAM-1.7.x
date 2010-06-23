/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "explicitSource.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(explicitSource, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        explicitSource,
        dictionary
    );
}

const Foam::wordList Foam::explicitSource::volumeModeTypeNames_
(
    IStringStream("(absolute specific)")()
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::explicitSource::setSelectedCellsFromPoints()
{
    labelHashSet selectedCells;

    forAll(points_, i)
    {
        label cellI = this->mesh().findCell(points_[i]);
        if (cellI >= 0)
        {
            selectedCells.insert(cellI);
        }

        label globalCellI = returnReduce(cellI, maxOp<label>());

        if (globalCellI < 0)
        {
            WarningIn("explicitSource::setSelectedCellsFromPoints()")
                << "Unable to find owner cell for point " << points_[i]
                << endl;
        }
    }

    this->cells() = selectedCells.toc();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::explicitSource::volumeModeType
Foam::explicitSource::wordToVolumeModeType
(
    const word& vmtName
) const
{
    forAll(volumeModeTypeNames_, i)
    {
        if (vmtName == volumeModeTypeNames_[i])
        {
            return volumeModeType(i);
        }
    }

    FatalErrorIn
    (
        "explicitSource<Type>::volumeModeType"
        "explicitSource<Type>::wordToVolumeModeType(const word&)"
    )   << "Unknown volumeMode type " << vmtName
        << ". Valid volumeMode types are:" << nl << volumeModeTypeNames_
        << exit(FatalError);

    return volumeModeType(0);
}


Foam::word Foam::explicitSource::volumeModeTypeToWord
(
    const volumeModeType& vmtType
) const
{
    if (vmtType > volumeModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return volumeModeTypeNames_[vmtType];
    }
}


void Foam::explicitSource::setFieldData(const dictionary& dict)
{
    scalarFields_.clear();
    vectorFields_.clear();

    wordList fieldTypes(dict.toc().size());
    wordList fieldNames(dict.toc().size());

    forAll(dict.toc(), i)
    {
        const word& fieldName = dict.toc()[i];
        IOobject io
        (
            fieldName,
            this->mesh().time().timeName(0),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );
        if (io.headerOk())
        {
            fieldTypes[i] = io.headerClassName();
            fieldNames[i] = dict.toc()[i];
        }
        else
        {
            FatalErrorIn
            (
                "explicitSource::setFieldData"
            )   << "header not OK " << io.name()
                << exit(FatalError);
        }
    }

    addField(scalarFields_, fieldTypes, fieldNames, dict);
    addField(vectorFields_, fieldTypes, fieldNames, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::explicitSource::explicitSource
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, dict, mesh),
    scalarFields_(0, *this),
    vectorFields_(0, *this),
    dict_(dict.subDict(typeName + "Coeffs")),
    volumeMode_(wordToVolumeModeType(dict_.lookup("volumeMode"))),
    points_(),
    volSource_(this->cells().size(), 1.0)
{
    setFieldData(dict_.subDict("fieldData"));

    // Set points if selectionMode is smPoints
    if (this->selectionMode() == smPoints)
    {
        dict_.lookup("points") >> points_;
        setSelectedCellsFromPoints();
        volSource_.setSize(points_.size(), 1.0);
    }

    const labelList& cellList = this->cells();
    scalar V = 0.0;
    if (volumeMode_ == vmAbsolute)
    {
        forAll(cellList, cellI)
        {
            volSource_[cellI] = mesh.V()[cellList[cellI]];
            V += volSource_[cellI];
        }
    }
    else
    {
        forAll(cellList, cellI)
        {
            V += mesh.V()[cellList[cellI]];
        }
    }

    reduce(V, sumOp<scalar>());

    Info<< "- selected " << returnReduce(cellList.size(), sumOp<label>())
        << " cell(s) with Volume: " << V << " in time activated sources "
        <<  endl;
}


void Foam::explicitSource::addSu(fvMatrix<scalar>& Eqn)
{
    Field<scalar>& source = Eqn.source();
    scalar data = scalarFields_[Eqn.psi().name()];
    addSources<scalar>(source, data);
}


void Foam::explicitSource::addSu(fvMatrix<vector>& Eqn)
{
    Field<vector>& source = Eqn.source();
    vector data = vectorFields_[Eqn.psi().name()];
    addSources<vector>(source, data);
}


void Foam::explicitSource::addSu(DimensionedField<scalar, volMesh>& field)
{
    scalar data = scalarFields_[field.name()];
    addSources<scalar>(field, data);
}


void Foam::explicitSource::addSu(DimensionedField<vector, volMesh>& field)
{
    vector data = vectorFields_[field.name()];
    addSources<vector>(field, data);
}


void Foam::explicitSource::addExplicitSources()
{
    scalarFields_.applySources();
    vectorFields_.applySources();
}


// ************************************************************************* //
