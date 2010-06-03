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


#include "exponential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable(pdf, exponential, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::exponential::exponential(const dictionary& dict, Random& rndGen)
:
    pdf(dict, rndGen),
    pdfDict_(dict.subDict(typeName + "PDF")),
    minValue_(readScalar(pdfDict_.lookup("minValue"))),
    maxValue_(readScalar(pdfDict_.lookup("maxValue"))),
    lambda_(pdfDict_.lookup("lambda")),
    ls_(lambda_),
    range_(maxValue_-minValue_)
{
    if (minValue_<0)
    {
        FatalErrorIn
        (
            "exponential::exponential(const dictionary& dict)"
        ) << " minValue = " << minValue_ << ", it must be >0." << abort(FatalError);
    }

    scalar sMax = 0;
    label n = lambda_.size();
    for (label i=0; i<n; i++)
    {
        scalar s = lambda_[i]*exp(-lambda_[i]*minValue_);
        for (label j=0; j<n; j++)
        {
            if (i!=j)
            {
                scalar y = lambda_[j]*exp(-lambda_[j]*minValue_);
                s += y;
            }
        }

        sMax = max(sMax, s);
    }

    for(label i=0; i<n; i++)
    {
        ls_[i] /= sMax;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::exponential::~exponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::exponential::sample() const
{
    scalar y = 0;
    scalar x = 0;
    label n = lambda_.size();
    bool success = false;

    while (!success)
    {
        x = minValue_ + range_*rndGen_.scalar01();
        y = rndGen_.scalar01();
        scalar p = 0.0;

        for(label i=0; i<n; i++)
        {
            p += ls_[i]*exp(-lambda_[i]*x);
        }

        if (y<p)
        {
            success = true;
        }
    }

    return x;
}


Foam::scalar Foam::exponential::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::exponential::maxValue() const
{
    return maxValue_;
}


// ************************************************************************* //
