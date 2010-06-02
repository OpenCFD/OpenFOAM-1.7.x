/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "topoSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(topoSet, 0);
defineRunTimeSelectionTable(topoSet, word);
defineRunTimeSelectionTable(topoSet, size);
defineRunTimeSelectionTable(topoSet, set);


// Construct named object from existing set.
autoPtr<topoSet> topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_
            ->find(setType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "topoSet::New(const word&, "
            "const polyMesh&, const word&, readOption, writeOption)"
        )   << "Unknown set type " << setType
            << endl << endl
            << "Valid set types : " << endl
            << wordConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<topoSet>(cstrIter()(mesh, name, r, w));
}


// Construct named object from size (non-existing set).
autoPtr<topoSet> topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
{
    sizeConstructorTable::iterator cstrIter =
        sizeConstructorTablePtr_
            ->find(setType);

    if (cstrIter == sizeConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "topoSet::New(const word&, "
            "const polyMesh&, const word&, const label, writeOption)"
        )   << "Unknown set type " << setType
            << endl << endl
            << "Valid set types : " << endl
            << sizeConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<topoSet>(cstrIter()(mesh, name, size, w));
}


// Construct named object from existing set.
autoPtr<topoSet> topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
{
    setConstructorTable::iterator cstrIter =
        setConstructorTablePtr_
            ->find(setType);

    if (cstrIter == setConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "topoSet::New(const word&, "
            "const polyMesh&, const word&, const topoSet&, writeOption)"
        )   << "Unknown set type " << setType
            << endl << endl
            << "Valid set types : " << endl
            << setConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<topoSet>(cstrIter()(mesh, name, set, w));
}


Foam::fileName topoSet::topoSet::localPath
(
    const polyMesh& mesh,
    const word& name
)
{
    return mesh.pointsInstance()/polyMesh::meshSubDir/"sets"/name;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update stored cell numbers using map.
// Do in two passes to prevent allocation if nothing changed.
void topoSet::topoSet::updateLabels(const labelList& map)
{
    // Iterate over map to see if anything changed
    bool changed = false;

    for
    (
        labelHashSet::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        if ((iter.key() < 0) || (iter.key() > map.size()))
        {
            FatalErrorIn
            (
                "topoSet::updateLabels(const labelList&, labelHashSet)"
            )   << "Illegal content " << iter.key() << " of set:" << name()
                << " of type " << type() << endl
                << "Value should be between 0 and " << map.size()-1
                << abort(FatalError);
        }

        label newCellI = map[iter.key()];

        if (newCellI != iter.key())
        {
            changed = true;

            break;
        }
    }

    // Relabel (use second Map to prevent overlapping)
    if (changed)
    {
        labelHashSet newSet(2*size());

        for
        (
            labelHashSet::const_iterator iter = begin();
            iter != end();
            ++iter
        )
        {
            label newCellI = map[iter.key()];

            if (newCellI >= 0)
            {
                newSet.insert(newCellI);
            }
        }

        transfer(newSet);
    }
}


void topoSet::topoSet::check(const label maxLabel)
{
    for
    (
        topoSet::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        if ((iter.key() < 0) || (iter.key() > maxLabel))
        {
            FatalErrorIn("topoSet::check(const label)")
                << "Illegal content " << iter.key() << " of set:" << name()
                << " of type " << type() << endl
                << "Value should be between 0 and " << maxLabel
                << abort(FatalError);
        }
    }
}


// Write maxElem elements, starting at iter. Updates iter and elemI.
void topoSet::writeDebug
(
    Ostream& os,
    const label maxElem,
    topoSet::const_iterator& iter,
    label& elemI
) const
{
    label n = 0;

    for (; (iter != end()) && (n < maxElem); ++iter)
    {
        if ((n != 0) && ((n % 10) == 0))
        {
            os << endl;
        }
        os << iter.key() << ' ';

        n++;
        elemI++;
    }
}


// Write maxElem elements, starting at iter. Updates iter and elemI.
void topoSet::writeDebug
(
    Ostream& os,
    const pointField& coords,
    const label maxElem,
    topoSet::const_iterator& iter,
    label& elemI
) const
{
    label n = 0;

    for (; (iter != end()) && (n < maxElem); ++iter)
    {
        if ((n != 0) && ((n % 3) == 0))
        {
            os << endl;
        }
        os << iter.key() << coords[iter.key()] << ' ';

        n++;
        elemI++;
    }
}


void topoSet::writeDebug
(
    Ostream& os,
    const pointField& coords,
    const label maxLen
) const
{
    // Bounding box of contents.
    boundBox bb(pointField(coords, toc()), true);

    os  << "Set bounding box: min = "
        << bb.min() << "    max = " << bb.max() << " meters. " << endl << endl;

    label n = 0;

    topoSet::const_iterator iter = begin();

    if (size() <= maxLen)
    {
        writeDebug(os, coords, maxLen, iter, n);
    }
    else
    {
        label halfLen = maxLen/2;

        os  << "Size larger than " << maxLen << ". Printing first and last "
            << halfLen << " elements:" << endl << endl;

        writeDebug(os, coords, halfLen, iter, n);

        os<< endl
          << "  .." << endl
          << endl;

        for (; n < size() - halfLen; ++n)
        {
            ++iter;
        }

        writeDebug(os, coords, halfLen, iter, n);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

topoSet::topoSet(const IOobject& obj, const word& wantedType)
:
    regIOobject(obj)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || (
            readOpt() == IOobject::READ_IF_PRESENT
         && headerOk()
        )
    )
    {
        if (readStream(wantedType).good())
        {
            readStream(wantedType) >> static_cast<labelHashSet&>(*this);

            close();
        }
    }
}


topoSet::topoSet
(
    const polyMesh& mesh,
    const word& wantedType,
    const word& name,
    readOption r,
    writeOption w
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.pointsInstance(),
            polyMesh::meshSubDir/"sets",
            mesh,
            r,
            w
        )
    )
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || (
            readOpt() == IOobject::READ_IF_PRESENT
         && headerOk()
        )
    )
    {
        if (readStream(wantedType).good())
        {
            readStream(wantedType) >> static_cast<labelHashSet&>(*this);

            close();
        }
    }
}


topoSet::topoSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.pointsInstance(),
            polyMesh::meshSubDir/"sets",
            mesh,
            NO_READ,
            w
        )
    ),
    labelHashSet(size)
{}


topoSet::topoSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& set,
    writeOption w
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.pointsInstance(),
            polyMesh::meshSubDir/"sets",
            mesh,
            NO_READ,
            w
        )
    ),
    labelHashSet(set)
{}


topoSet::topoSet(const IOobject& obj, const label size)
:
    regIOobject(obj),
    labelHashSet(size)
{}


topoSet::topoSet(const IOobject& obj, const labelHashSet& set)
:
    regIOobject(obj),
    labelHashSet(set)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

topoSet::~topoSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void topoSet::invert(const label maxLen)
{
    // Keep copy of current set.
    labelHashSet currentSet(*this);

    clear();
    resize(2*(maxLen - currentSet.size()));

    for (label cellI = 0; cellI < maxLen; cellI++)
    {
        if (!currentSet.found(cellI))
        {
            insert(cellI);
        }
    }

}


void topoSet::subset(const topoSet& set)
{
    // Keep copy of current set.
    labelHashSet currentSet(*this);

    clear();
    resize(2*min(currentSet.size(), set.size()));

    for
    (
        labelHashSet::const_iterator iter = currentSet.begin();
        iter != currentSet.end();
        ++iter
    )
    {
        if (set.found(iter.key()))
        {
            // element present in both currentSet and set.
            insert(iter.key());
        }
    }
}


void topoSet::addSet(const topoSet& set)
{
    for
    (
        topoSet::const_iterator iter = set.begin();
        iter != set.end();
        ++iter
    )
    {
        insert(iter.key());
    }
}


void topoSet::deleteSet(const topoSet& set)
{
    for
    (
        topoSet::const_iterator iter = set.begin();
        iter != set.end();
        ++iter
    )
    {
        erase(iter.key());
    }
}


void topoSet::sync(const polyMesh&)
{
    notImplemented("topoSet::sync(const polyMesh&)");
}


void topoSet::writeDebug(Ostream& os, const label maxLen) const
{
    label n = 0;

    topoSet::const_iterator iter = begin();

    if (size() <= maxLen)
    {
        writeDebug(os, maxLen, iter, n);
    }
    else
    {
        label halfLen = maxLen/2;

        os  << "Size larger than " << maxLen << ". Printing first and last "
            << halfLen << " elements:" << endl << endl;

        writeDebug(os, halfLen, iter, n);

        os<< endl
          << "  .." << endl
          << endl;

        for (; n < size() - halfLen; ++n)
        {
            ++iter;
        }

        writeDebug(os, halfLen, iter, n);
    }
}


void topoSet::writeDebug
(
    Ostream&,
    const primitiveMesh&,
    const label
) const
{
    notImplemented
    (
        "topoSet::writeDebug(Ostream&, const primitiveMesh&, const label)"
    );
}


bool topoSet::writeData(Ostream& os) const
{
    return (os << *this).good();
}


void topoSet::updateMesh(const mapPolyMesh&)
{
    notImplemented("topoSet::updateMesh(const mapPolyMesh&)");
}


//- Return max index+1.
label topoSet::maxSize(const polyMesh&) const
{
    notImplemented("topoSet::maxSize(const polyMesh&)");

    return -1;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void topoSet::operator=(const topoSet& rhs)
{
    labelHashSet::operator=(rhs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
