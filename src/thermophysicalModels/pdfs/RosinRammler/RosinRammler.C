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

#include "RosinRammler.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace pdfs
    {
        defineTypeNameAndDebug(RosinRammler, 0);
        addToRunTimeSelectionTable(pdf, RosinRammler, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdfs::RosinRammler::RosinRammler(const dictionary& dict, Random& rndGen)
:
    pdf(typeName, dict, rndGen),
    minValue_(readScalar(pdfDict_.lookup("minValue"))),
    maxValue_(readScalar(pdfDict_.lookup("maxValue"))),
    d_(readScalar(pdfDict_.lookup("d"))),
    n_(readScalar(pdfDict_.lookup("n")))
{
    check();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdfs::RosinRammler::~RosinRammler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::pdfs::RosinRammler::sample() const
{

    scalar K = 1.0 - exp(-pow((maxValue_ - minValue_)/d_, n_));
    scalar y = rndGen_.scalar01();
    scalar x = minValue_ + d_*::pow(-log(1.0 - y*K), 1.0/n_);
    return x;
}


Foam::scalar Foam::pdfs::RosinRammler::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::pdfs::RosinRammler::maxValue() const
{
    return maxValue_;
}


// ************************************************************************* //
