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

Description
    A subset of mesh points.

\*---------------------------------------------------------------------------*/

#include "pointZone.H"
#include "addToRunTimeSelectionTable.H"
#include "pointZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointZone, 0);
    defineRunTimeSelectionTable(pointZone, dictionary);
    addToRunTimeSelectionTable(pointZone, pointZone, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Map<Foam::label>& Foam::pointZone::pointLookupMap() const
{
    if (!pointLookupMapPtr_)
    {
        calcPointLookupMap();
    }

    return *pointLookupMapPtr_;
}


void Foam::pointZone::calcPointLookupMap() const
{
    if (debug)
    {
        Info<< "void pointZone::calcPointLookupMap() const : "
            << "Calculating point lookup map"
            << endl;
    }

    if (pointLookupMapPtr_)
    {
        FatalErrorIn
        (
            "void pointZone::calcPointLookupMap() const"
        )   << "point lookup map already calculated"
            << abort(FatalError);
    }

    const labelList& addr = *this;

    pointLookupMapPtr_ = new Map<label>(2*addr.size());
    Map<label>& plm = *pointLookupMapPtr_;

    forAll (addr, pointI)
    {
        plm.insert(addr[pointI], pointI);
    }

    if (debug)
    {
        Info<< "void pointZone::calcPointLookupMap() const : "
            << "Finished calculating point lookup map"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pointZone::pointZone
(
    const word& name,
    const labelList& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    labelList(addr),
    name_(name),
    index_(index),
    zoneMesh_(zm),
    pointLookupMapPtr_(NULL)
{}


Foam::pointZone::pointZone
(
    const word& name,
    const Xfer<labelList>& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    labelList(addr),
    name_(name),
    index_(index),
    zoneMesh_(zm),
    pointLookupMapPtr_(NULL)
{}


// Construct from dictionary
Foam::pointZone::pointZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const pointZoneMesh& zm
)
:
    labelList(dict.lookup("pointLabels")),
    name_(name),
    index_(index),
    zoneMesh_(zm),
    pointLookupMapPtr_(NULL)
{}


// Construct given the original zone and resetting the
// point list and zone mesh information
Foam::pointZone::pointZone
(
    const pointZone& pz,
    const labelList& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    labelList(addr),
    name_(pz.name()),
    index_(index),
    zoneMesh_(zm),
    pointLookupMapPtr_(NULL)
{}


Foam::pointZone::pointZone
(
    const pointZone& pz,
    const Xfer<labelList>& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    labelList(addr),
    name_(pz.name()),
    index_(index),
    zoneMesh_(zm),
    pointLookupMapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointZone::~pointZone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::pointZone::whichPoint(const label globalPointID) const
{
    const Map<label>& plm = pointLookupMap();

    Map<label>::const_iterator plmIter = plm.find(globalPointID);

    if (plmIter == plm.end())
    {
        return -1;
    }
    else
    {
        return plmIter();
    }
}


const Foam::pointZoneMesh& Foam::pointZone::zoneMesh() const
{
    return zoneMesh_;
}


void Foam::pointZone::clearAddressing()
{
    deleteDemandDrivenData(pointLookupMapPtr_);
}


bool Foam::pointZone::checkDefinition(const bool report) const
{
    const labelList& addr = *this;

    bool boundaryError = false;

    forAll(addr, i)
    {
        if (addr[i] < 0 || addr[i] >= zoneMesh_.mesh().points().size())
        {
            boundaryError = true;

            if (report)
            {
                SeriousErrorIn
                (
                    "bool pointZone::checkDefinition("
                    "const bool report) const"
                )   << "Zone " << name()
                    << " contains invalid point label " << addr[i] << nl
                    << "Valid point labels are 0.."
                    << zoneMesh_.mesh().points().size()-1 << endl;
            }
        }
    }
    return boundaryError;
}


void Foam::pointZone::write(Ostream& os) const
{
    os  << nl << name()
        << nl << static_cast<const labelList&>(*this);
}


void Foam::pointZone::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry("pointLabels", os);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pointZone::operator=(const pointZone& cz)
{
    clearAddressing();
    labelList::operator=(cz);
}


void Foam::pointZone::operator=(const labelList& addr)
{
    clearAddressing();
    labelList::operator=(addr);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pointZone& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const pointZone& p");
    return os;
}


// ************************************************************************* //
