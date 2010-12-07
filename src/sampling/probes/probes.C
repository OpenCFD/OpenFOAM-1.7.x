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

#include "probes.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(probes, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::probes::findElements(const fvMesh& mesh)
{
    if (elementList_.empty())
    {
        elementList_.setSize(probeLocations_.size());

        forAll(probeLocations_, probeI)
        {
            elementList_[probeI] = mesh.findCell(probeLocations_[probeI]);

            if (debug && elementList_[probeI] != -1)
            {
                Pout<< "probes : found point " << probeLocations_[probeI]
                    << " in cell " << elementList_[probeI] << endl;
            }
        }


        // Check if all probes have been found.
        forAll(elementList_, probeI)
        {
            label cellI = elementList_[probeI];

            // Check at least one processor with cell.
            reduce(cellI, maxOp<label>());

            if (cellI == -1)
            {
                if (Pstream::master())
                {
                    WarningIn("probes::read()")
                        << "Did not find location " << probeLocations_[probeI]
                        << " in any cell. Skipping location." << endl;
                }
            }
            else
            {
                // Make sure location not on two domains.
                if (elementList_[probeI] != -1 && elementList_[probeI] != cellI)
                {
                    WarningIn("probes::read()")
                        << "Location " << probeLocations_[probeI]
                        << " seems to be on multiple domains:"
                        << " cell " << elementList_[probeI]
                        << " on my domain " << Pstream::myProcNo()
                        << " and cell " << cellI << " on some other domain."
                        << endl
                        << "This might happen if the probe location is on"
                        << " a processor patch. Change the location slightly"
                        << " to prevent this." << endl;
                }
            }
        }
    }
}


bool Foam::probes::checkFieldTypes()
{
    wordList fieldTypes(fieldNames_.size());

    // check files for a particular time
    if (loadFromFiles_)
    {
        forAll(fieldNames_, fieldI)
        {
            IOobject io
            (
                fieldNames_[fieldI],
                obr_.time().timeName(),
                refCast<const polyMesh>(obr_),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (io.headerOk())
            {
                fieldTypes[fieldI] = io.headerClassName();
            }
            else
            {
                fieldTypes[fieldI] = "(notFound)";
            }
        }
    }
    else
    {
        // check objectRegistry
        forAll(fieldNames_, fieldI)
        {
            objectRegistry::const_iterator iter =
                obr_.find(fieldNames_[fieldI]);

            if (iter != obr_.end())
            {
                fieldTypes[fieldI] = iter()->type();
            }
            else
            {
                fieldTypes[fieldI] = "(notFound)";
            }
        }
    }


    label nFields = 0;

    // classify fieldTypes
    nFields += countFields(scalarFields_, fieldTypes);
    nFields += countFields(vectorFields_, fieldTypes);
    nFields += countFields(sphericalTensorFields_, fieldTypes);
    nFields += countFields(symmTensorFields_, fieldTypes);
    nFields += countFields(tensorFields_, fieldTypes);

    // concatenate all the lists into foundFields
    wordList foundFields(nFields);

    label fieldI = 0;
    forAll(scalarFields_, i)
    {
        foundFields[fieldI++] = scalarFields_[i];
    }
    forAll(vectorFields_, i)
    {
        foundFields[fieldI++] = vectorFields_[i];
    }
    forAll(sphericalTensorFields_, i)
    {
        foundFields[fieldI++] = sphericalTensorFields_[i];
    }
    forAll(symmTensorFields_, i)
    {
        foundFields[fieldI++] = symmTensorFields_[i];
    }
    forAll(tensorFields_, i)
    {
        foundFields[fieldI++] = tensorFields_[i];
    }

    if (Pstream::master())
    {
        fileName probeDir;

        fileName probeSubDir = name_;

        if (obr_.name() != polyMesh::defaultRegion)
        {
            probeSubDir = probeSubDir/obr_.name();
        }
        probeSubDir = probeSubDir/obr_.time().timeName();

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            probeDir = obr_.time().path()/".."/probeSubDir;
        }
        else
        {
            probeDir = obr_.time().path()/probeSubDir;
        }

        // Close the file if any fields have been removed.
        forAllIter(HashPtrTable<OFstream>, probeFilePtrs_, iter)
        {
            if (findIndex(foundFields, iter.key()) == -1)
            {
                if (debug)
                {
                    Pout<< "close stream: " << iter()->name() << endl;
                }

                delete probeFilePtrs_.remove(iter);
            }
        }

        // Open new files for new fields. Keep existing files.

        probeFilePtrs_.resize(2*foundFields.size());

        forAll(foundFields, fieldI)
        {
            const word& fldName = foundFields[fieldI];

            // Check if added field. If so open a stream for it.

            if (!probeFilePtrs_.found(fldName))
            {
                // Create directory if does not exist.
                mkDir(probeDir);

                OFstream* sPtr = new OFstream(probeDir/fldName);

                if (debug)
                {
                    Pout<< "open  stream: " << sPtr->name() << endl;
                }

                probeFilePtrs_.insert(fldName, sPtr);

                unsigned int w = IOstream::defaultPrecision() + 7;

                for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                {
                    *sPtr<< '#' << setw(IOstream::defaultPrecision() + 6)
                        << vector::componentNames[cmpt];

                    forAll(probeLocations_, probeI)
                    {
                        *sPtr<< ' ' << setw(w) << probeLocations_[probeI][cmpt];
                    }
                    *sPtr << endl;
                }

                *sPtr<< '#' << setw(IOstream::defaultPrecision() + 6)
                    << "Time" << endl;
            }
        }

        if (debug)
        {
            Pout<< "Probing fields:" << foundFields << nl
                << "Probing locations:" << probeLocations_ << nl
                << endl;
        }
    }


    return nFields > 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::probes::probes
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    loadFromFiles_(loadFromFiles),
    fieldNames_(0),
    probeLocations_(0),
    scalarFields_(),
    vectorFields_(),
    sphericalTensorFields_(),
    symmTensorFields_(),
    tensorFields_(),
    elementList_(0),
    probeFilePtrs_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::probes::~probes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::probes::execute()
{
    // Do nothing - only valid on write
}


void Foam::probes::end()
{
    // Do nothing - only valid on write
}


void Foam::probes::write()
{
    if (probeLocations_.size() && checkFieldTypes())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);
    }
}


void Foam::probes::read(const dictionary& dict)
{
    dict.lookup("fields") >> fieldNames_;
    dict.lookup("probeLocations") >> probeLocations_;

    // Force all cell locations to be redetermined
    elementList_.clear();
    findElements(refCast<const fvMesh>(obr_));
    checkFieldTypes();
}


// ************************************************************************* //
