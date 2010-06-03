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

\*----------------------------------------------------------------------------*/

#include "ensightPart.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
   defineTypeNameAndDebug(ensightPart, 0);
   defineTemplateTypeNameAndDebug(IOPtrList<ensightPart>, 0);
   defineRunTimeSelectionTable(ensightPart, istream);
}

Foam::List<Foam::word> Foam::ensightPart::elemTypes_(0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::ensightPart::isFieldDefined(const List<scalar>& field) const
{
    forAll(elemLists_, elemI)
    {
        const labelList& idList = elemLists_[elemI];

        forAll(idList, i)
        {
            label id = idList[i];

            if (id >= field.size() || isnan(field[id]))
            {
                return false;
            }
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPart::ensightPart
()
:
    number_(0),
    name_(""),
    elemLists_(0),
    offset_(0),
    size_(0),
    isCellData_(true),
    matId_(0),
    meshPtr_(0)
{}


Foam::ensightPart::ensightPart
(
    label partNumber,
    const string& partDescription
)
:
    number_(partNumber),
    name_(partDescription),
    elemLists_(0),
    offset_(0),
    size_(0),
    isCellData_(true),
    matId_(0),
    meshPtr_(0)
{}


Foam::ensightPart::ensightPart
(
    label partNumber,
    const string& partDescription,
    const polyMesh& pMesh
)
:
    number_(partNumber),
    name_(partDescription),
    elemLists_(0),
    offset_(0),
    size_(0),
    isCellData_(true),
    matId_(0),
    meshPtr_(&pMesh)
{}


Foam::ensightPart::ensightPart(const ensightPart& part)
:
    number_(part.number_),
    name_(part.name_),
    elemLists_(part.elemLists_),
    offset_(part.offset_),
    size_(part.size_),
    isCellData_(part.isCellData_),
    matId_(part.matId_),
    meshPtr_(part.meshPtr_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ensightPart> Foam::ensightPart::New(Istream& is)
{
    word partType(is);

    istreamConstructorTable::iterator cstrIter =
        istreamConstructorTablePtr_->find(partType);

    if (cstrIter == istreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "ensightPart::New(Istream&)",
            is
        )   << "unknown ensightPart type " << partType << endl << endl
            << "Valid ensightPart types are :" << endl
            << istreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<ensightPart>(cstrIter()(is));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPart::~ensightPart()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightPart::reconstruct(Istream& is)
{
    dictionary dict(is);
    dict.lookup("id") >> number_;
    dict.lookup("name") >> name_;
    dict.readIfPresent("offset", offset_);

    // populate elemLists_
    elemLists_.setSize(elementTypes().size());

    forAll(elementTypes(), elemI)
    {
        word key(elementTypes()[elemI]);

        elemLists_[elemI].clear();
        dict.readIfPresent(key, elemLists_[elemI]);

        size_ += elemLists_[elemI].size();
    }

    is.check("ensightPart::reconstruct(Istream&)");
}


void Foam::ensightPart::renumber(labelList const& origId)
{
    // transform to global values first
    if (offset_)
    {
        forAll(elemLists_, elemI)
        {
            labelList& idList = elemLists_[elemI];
            forAll(idList, i)
            {
                idList[i] += offset_;
            }
        }

        offset_ = 0;
    }

    if (origId.size())
    {
        forAll(elemLists_, elemI)
        {
            inplaceRenumber(origId, elemLists_[elemI]);
        }
    }
}


// ************************************************************************* //
