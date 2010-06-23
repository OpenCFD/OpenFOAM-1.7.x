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

#include "basicSource.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicSource, 0);
    defineRunTimeSelectionTable(basicSource, dictionary);
}


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

const Foam::wordList Foam::basicSource::selectionModeTypeNames_
(
    IStringStream("(points cellSet cellZone all)")()
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::basicSource::selectionModeType Foam::basicSource::wordToSelectionModeType
(
    const word& smtName
) const
{
    forAll(selectionModeTypeNames_, i)
    {
        if (smtName == selectionModeTypeNames_[i])
        {
            return selectionModeType(i);
        }
    }

    FatalErrorIn
    (
        "basicSource::selectionModeType"
        "basicSource::wordToSelectionModeType"
        "("
            "const word&"
        ")"
    )   << "Unknown selectionMode type " << smtName
        << ". Valid selectionMode types are:" << nl << selectionModeTypeNames_
        << exit(FatalError);

    return selectionModeType(0);
}


Foam::word Foam::basicSource::selectionModeTypeToWord
(
    const selectionModeType& smtType
) const
{
    if (smtType > selectionModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return selectionModeTypeNames_[smtType];
    }
}


void Foam::basicSource::setSelection(const dictionary& dict)
{
    switch (selectionMode_)
    {
        case smPoints:
        {
            // Do nothing. It should be sorted out by derived class
            break;
        }
        case smCellSet:
        {
            dict.lookup("cellSet") >> cellSetName_;
            break;
        }
        case smCellZone:
        {
            dict.lookup("cellZone") >> cellSetName_;
            break;
        }
        case smAll:
        {
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "basicSource::setSelection(const dictionary&)"
            )   << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }
}


void Foam::basicSource::setCellSet()
{
    Info<< incrIndent << indent << "Source: " << name_ << endl;
    switch (selectionMode_)
    {
        case smPoints:
        {
            break;
        }
        case smCellSet:
        {
            Info<< indent << "- selecting cells using cellSet "
                << cellSetName_ << endl;

            cellSet selectedCells(mesh_, cellSetName_);
            cells_ = selectedCells.toc();

            break;
        }
        case smCellZone:
        {
            Info<< indent << "- selecting cells using cellZone "
                << cellSetName_ << endl;
            label zoneID = mesh_.cellZones().findZoneID(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorIn("basicSource<Type>::setCellIds()")
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh_.cellZones().names()
                    << exit(FatalError);
            }
            cells_ = mesh_.cellZones()[zoneID];

            break;
        }
        case smAll:
        {
            Info<< indent << "- selecting all cells" << endl;
            cells_ = identity(mesh_.nCells());

            break;
        }
        default:
        {
            FatalErrorIn("basicSource<Type>::setCellIds()")
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    // Set volume information
    if (selectionMode_ != smPoints)
    {
        V_ = 0.0;
        forAll(cells_, i)
        {
            V_ += mesh_.V()[cells_[i]];
        }
        reduce(V_, sumOp<scalar>());

        Info<< indent << "- selected "
            << returnReduce(cells_.size(), sumOp<label>())
            << " cell(s) with volume " << V_ << nl << decrIndent << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSource::basicSource
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    active_(readBool(dict_.lookup("active"))),
    timeStart_(readScalar(dict_.lookup("timeStart"))),
    duration_(readScalar(dict_.lookup("duration"))),
    selectionMode_(wordToSelectionModeType(dict_.lookup("selectionMode"))),
    cellSetName_("none"),
    V_(1.0)
{
    setSelection(dict_);

    setCellSet();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicSource> Foam::basicSource::New
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word typeModel(dict.lookup("typeModel"));

    Info<< "Selecting model type " << typeModel << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(typeModel);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicSource::New(const volVectorField&, "
            "const surfaceScalarField&, transportModel&)"
        )   << "Unknown Model type " << typeModel
            << nl << nl
            << "Valid model types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<basicSource>(cstrIter()(name, dict, mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::basicSource::isActive()
{
    if
    (
        active_
     && (mesh_.time().value() >= timeStart_)
     && (mesh_.time().value() <= timeEnd())
    )
    {
        // Update the cell set if the mesh is changing
        if (mesh_.changing())
        {
            setCellSet();
        }
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
