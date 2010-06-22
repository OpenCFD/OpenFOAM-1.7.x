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

#include "fixedValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace pdfs
    {
        defineTypeNameAndDebug(fixedValue, 0);
        addToRunTimeSelectionTable(pdf, fixedValue, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdfs::fixedValue::fixedValue(const dictionary& dict, Random& rndGen)
:
    pdf(typeName, dict, rndGen),
    value_(readScalar(pdfDict_.lookup("value")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdfs::fixedValue::~fixedValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::pdfs::fixedValue::fixedValue::sample() const
{
    return value_;
}


Foam::scalar Foam::pdfs::fixedValue::fixedValue::minValue() const
{
    return -VGREAT;
}


Foam::scalar Foam::pdfs::fixedValue::fixedValue::maxValue() const
{
    return VGREAT;
}


// ************************************************************************* //
