/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "fieldValue.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldValue, 0);

    defineTemplateTypeNameAndDebug(IOList<vector>, 0);
    defineTemplateTypeNameAndDebug(IOList<sphericalTensor>, 0);
    defineTemplateTypeNameAndDebug(IOList<symmTensor>, 0);
    defineTemplateTypeNameAndDebug(IOList<tensor>, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValue::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldValue::movePoints(const Field<point>&)
{
    // Do nothing
}


void Foam::fieldValue::makeFile()
{
    // Create the forces file if not already created
    if (outputFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating output file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName outputDir;
            word startTimeName =
                obr_.time().timeName(obr_.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                outputDir =
                    obr_.time().path()/".."/name_/startTimeName;
            }
            else
            {
                outputDir = obr_.time().path()/name_/startTimeName;
            }

            // Create directory if does not exist
            mkDir(outputDir);

            // Open new file at start up
            outputFilePtr_.reset(new OFstream(outputDir/(type() + ".dat")));

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::fieldValue::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);
        dict.lookup("fields") >> fields_;
        dict.lookup("valueOutput") >> valueOutput_;
    }
}


void Foam::fieldValue::write()
{
    if (active_)
    {
        if (log_)
        {
            Info<< type() << " " << name_ << " output:" << nl;
        }

        makeFile();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValue::fieldValue
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
    sourceName_(dict.lookup("sourceName")),
    fields_(dict.lookup("fields")),
    valueOutput_(dict.lookup("valueOutput")),
    outputFilePtr_(NULL)
{
    // Only active if obr is an fvMesh
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        WarningIn
        (
            "fieldValue::fieldValue"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << nl << endl;
        active_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValue::execute()
{
    // Do nothing
}


void Foam::fieldValue::end()
{
    // Do nothing
}


// ************************************************************************* //
