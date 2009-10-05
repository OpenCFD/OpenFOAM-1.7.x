/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "absorptionCoeffs.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::absorptionCoeffs::absorptionCoeffs(Istream& is)
:
   Tcommon_(readScalar(is)),
   Tlow_(readScalar(is)),
   Thigh_(readScalar(is)),
   invTemp_(readBool(is))
{
    for (label coefLabel=0; absorptionCoeffs::nCoeffs_; coefLabel++)
    {
        is >> highACoeffs_[coefLabel];
    }

    for (label coefLabel=0; absorptionCoeffs::nCoeffs_; coefLabel++)
    {
        is >> lowACoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::radiation::absorptionCoeffs::~absorptionCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::absorptionCoeffs::checkT(const scalar T) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        FatalErrorIn
        (
            "absorptionCoeffs::checkT(const scalar T) const"
        )   << "attempt to use absCoeff out of temperature range:" << nl
            << "    " << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << nl << abort(FatalError);
    }
}


const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::coeffs
(
    const scalar T
) const
{
    checkT(T);

    if (T < Tcommon_)
    {
        return lowACoeffs_;
    }
    else
    {
        return highACoeffs_;
    }
}


void Foam::radiation::absorptionCoeffs::initialise(Istream&)
{
    absorptionCoeffs(Istream);
}


void Foam::radiation::absorptionCoeffs::initialise(const dictionary& dict)
{
    dict.lookup("Tcommon") >> Tcommon_;
    dict.lookup("Tlow") >> Tlow_;
    dict.lookup("Thigh") >> Thigh_;
    dict.lookup("invTemp") >> invTemp_;

    dict.lookup("loTcoeffs") >> lowACoeffs_;
    dict.lookup("hiTcoeffs") >> highACoeffs_;
}


// ************************************************************************* //
