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

#include "pdf.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace pdfs
    {
        defineTypeNameAndDebug(pdf, 0);
        defineRunTimeSelectionTable(pdf, dictionary);
    }
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::pdfs::pdf::check() const
{
    if (minValue() < 0)
    {
        FatalErrorIn("pdfs::pdf::check() const")
            << type() << "PDF: Minimum value must be greater than zero." << nl
            << "Supplied minValue = " << minValue()
            << abort(FatalError);
    }

    if (maxValue() < minValue())
    {
        FatalErrorIn("pdfs::pdf::check() const")
            << type() << "PDF: Maximum value is smaller than the minimum value:"
            << nl << "    maxValue = " << maxValue()
            << ", minValue = " << minValue()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdfs::pdf::pdf(const word& name, const dictionary& dict, Random& rndGen)
:
    pdfDict_(dict.subDict(name + "PDF")),
    rndGen_(rndGen)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdfs::pdf::~pdf()
{}


// ************************************************************************* //
