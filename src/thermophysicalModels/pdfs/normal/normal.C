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

#include "normal.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace pdfs
    {
        defineTypeNameAndDebug(normal, 0);
        addToRunTimeSelectionTable(pdf, normal, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdfs::normal::normal(const dictionary& dict, Random& rndGen)
:
    pdf(typeName, dict, rndGen),
    minValue_(readScalar(pdfDict_.lookup("minValue"))),
    maxValue_(readScalar(pdfDict_.lookup("maxValue"))),
    expectation_(readScalar(pdfDict_.lookup("expectation"))),
    variance_(readScalar(pdfDict_.lookup("variance"))),
    a_(0.147)
{
    if (minValue_ < 0)
    {
        FatalErrorIn("normal::normal(const dictionary&, Random&)")
            << "Minimum value must be greater than zero. "
            << "Supplied minValue = " << minValue_
            << abort(FatalError);
    }

    if (maxValue_ < minValue_)
    {
        FatalErrorIn("normal::normal(const dictionary&, Random&)")
            << "Maximum value is smaller than the minimum value:"
            << "    maxValue = " << maxValue_ << ", minValue = " << minValue_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdfs::normal::~normal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::pdfs::normal::sample() const
{

    scalar a = erf((minValue_ - expectation_)/variance_);
    scalar b = erf((maxValue_ - expectation_)/variance_);

    scalar y = rndGen_.scalar01();
    scalar x = erfInv(y*(b - a) + a)*variance_ + expectation_;

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    x = min(max(x, minValue_), maxValue_);

    return x;
}


Foam::scalar Foam::pdfs::normal::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::pdfs::normal::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::pdfs::normal::erfInv(const scalar y) const
{
    scalar k = 2.0/(mathematicalConstant::pi*a_) +  0.5*log(1.0 - y*y);
    scalar h = log(1.0 - y*y)/a_;
    scalar x = sqrt(-k + sqrt(k*k - h));
    if (y < 0.0)
    {
        x *= -1.0;
    }
    return x;
}


// ************************************************************************* //
