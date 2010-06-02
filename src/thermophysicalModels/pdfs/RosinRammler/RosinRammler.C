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

#include "RosinRammler.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RosinRammler, 0);

    addToRunTimeSelectionTable(pdf, RosinRammler, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RosinRammler::RosinRammler(const dictionary& dict, Random& rndGen)
:
    pdf(dict, rndGen),
    pdfDict_(dict.subDict(typeName + "PDF")),
    minValue_(readScalar(pdfDict_.lookup("minValue"))),
    maxValue_(readScalar(pdfDict_.lookup("maxValue"))),
    d_(pdfDict_.lookup("d")),
    n_(pdfDict_.lookup("n")),
    ls_(d_),
    range_(maxValue_-minValue_)
{
    if (minValue_<0)
    {
        FatalErrorIn
        (
            "RosinRammler::RosinRammler(const dictionary& dict)"
        ) << " minValue = " << minValue_ << ", it must be >0." << abort(FatalError);
    }

    if (maxValue_<minValue_)
    {
        FatalErrorIn
        (
            "RosinRammler::RosinRammler(const dictionary& dict)"
        ) << " maxValue is smaller than minValue." << abort(FatalError);
    }

    // find max value so that it can be normalized to 1.0
    scalar sMax = 0;
    label n = d_.size();
    for (label i=0; i<n; i++)
    {
        scalar s = exp(-1.0);
        for (label j=0; j<n; j++)
        {
            if (i!=j)
            {
                scalar xx = pow(d_[j]/d_[i], n_[j]);
                scalar y = xx*exp(-xx);
                s += y;
            }
        }

        sMax = max(sMax, s);
    }

    for (label i=0; i<n; i++)
    {
        ls_[i] /= sMax;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RosinRammler::~RosinRammler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::RosinRammler::sample() const
{
    scalar y = 0;
    scalar x = 0;
    label n = d_.size();
    bool success = false;

    while (!success)
    {
        x = minValue_ + range_*rndGen_.scalar01();
        y = rndGen_.scalar01();
        scalar p = 0.0;

        for (label i=0; i<n; i++)
        {
            scalar xx = pow(x/d_[i], n_[i]);
            p += ls_[i]*xx*exp(-xx);
        }

        if (y<p)
        {
            success = true;
        }
    }

    return x;
}


Foam::scalar Foam::RosinRammler::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::RosinRammler::maxValue() const
{
    return maxValue_;
}


// ************************************************************************* //
