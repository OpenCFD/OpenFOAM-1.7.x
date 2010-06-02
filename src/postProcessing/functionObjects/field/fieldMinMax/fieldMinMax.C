/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "fieldMinMax.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldMinMax, 0);
}


template<>
const char* Foam::NamedEnum<Foam::fieldMinMax::modeType, 2>::names[] =
{
    "magnitude",
    "component"
};


const Foam::NamedEnum<Foam::fieldMinMax::modeType, 2>
Foam::fieldMinMax::modeTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldMinMax::fieldMinMax
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    log_(false),
    mode_(mdMag),
    fieldSet_(),
    fieldMinMaxFilePtr_(NULL)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "fieldMinMax::fieldMinMax"
            "(const objectRegistry& obr, const dictionary& dict)"
        )   << "No fvMesh available, deactivating."
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldMinMax::~fieldMinMax()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldMinMax::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);

        mode_ = modeTypeNames_[dict.lookup("mode")];
        dict.lookup("fields") >> fieldSet_;
    }
}


void Foam::fieldMinMax::makeFile()
{
    // Create the fieldMinMax file if not already created
    if (fieldMinMaxFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating fieldMinMax file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName fieldMinMaxDir;
            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                fieldMinMaxDir =
                    obr_.time().path()/".."/name_/obr_.time().timeName();
            }
            else
            {
                fieldMinMaxDir =
                    obr_.time().path()/name_/obr_.time().timeName();
            }

            // Create directory if does not exist.
            mkDir(fieldMinMaxDir);

            // Open new file at start up
            fieldMinMaxFilePtr_.reset
            (
                new OFstream(fieldMinMaxDir/(type() + ".dat"))
            );

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::fieldMinMax::writeFileHeader()
{
    if (fieldMinMaxFilePtr_.valid())
    {
        fieldMinMaxFilePtr_()
            << "# Time" << tab << "field" << tab << "min" << tab << "max"
            << endl;
    }
}


void Foam::fieldMinMax::execute()
{
    // Do nothing - only valid on write
}


void Foam::fieldMinMax::end()
{
    // Do nothing - only valid on write
}


void Foam::fieldMinMax::write()
{
    if (active_)
    {
        // Create the fieldMinMax file if not already created
        makeFile();

        forAll(fieldSet_, fieldI)
        {
            calcMinMaxFields<scalar>(fieldSet_[fieldI]);
            calcMinMaxFields<vector>(fieldSet_[fieldI]);
            calcMinMaxFields<sphericalTensor>(fieldSet_[fieldI]);
            calcMinMaxFields<symmTensor>(fieldSet_[fieldI]);
            calcMinMaxFields<tensor>(fieldSet_[fieldI]);
        }
    }
}


template<>
void Foam::fieldMinMax::calcMinMaxFields<Foam::scalar>
(
    const word& fieldName
)
{
    if (obr_.foundObject<volScalarField>(fieldName))
    {
        const volScalarField& field =
            obr_.lookupObject<volScalarField>(fieldName);
        scalar minValue = min(field).value();
        scalar maxValue = max(field).value();

        if (Pstream::master())
        {
            fieldMinMaxFilePtr_() << obr_.time().value() << tab
                << fieldName << tab << minValue << tab << maxValue << endl;

            if (log_)
            {
                Info<< "fieldMinMax output:" << nl
                    << "    min(" << fieldName << ") = " << minValue << nl
                    << "    max(" << fieldName << ") = " << maxValue << nl
                    << endl;
            }
        }
    }
}


// ************************************************************************* //
