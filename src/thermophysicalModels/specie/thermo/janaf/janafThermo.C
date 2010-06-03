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

#include "janafThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
Foam::janafThermo<equationOfState>::janafThermo(Istream& is)
:
    equationOfState(is),
    Tlow_(readScalar(is)),
    Thigh_(readScalar(is)),
    Tcommon_(readScalar(is))
{
    if (Tlow_ >= Thigh_)
    {
        FatalIOErrorIn
        (
            "janafThermo<equationOfState>::janafThermo(Istream& is)",
            is
        )   << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalIOError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalIOErrorIn
        (
            "janafThermo<equationOfState>::janafThermo(Istream& is)",
            is
        )   << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalIOError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalIOErrorIn
        (
            "janafThermo<equationOfState>::janafThermo(Istream& is)",
            is
        )   << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalIOError);
    }

    for
    (
        register label coefLabel=0;
        coefLabel<janafThermo<equationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        is >> highCpCoeffs_[coefLabel];
    }

    for
    (
        register label coefLabel=0;
        coefLabel<janafThermo<equationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        is >> lowCpCoeffs_[coefLabel];
    }

    // Check state of Istream
    is.check("janafThermo::janafThermo(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const janafThermo<equationOfState>& jt
)
{
    os  << static_cast<const equationOfState&>(jt) << nl
        << "    " << jt.Tlow_
        << tab << jt.Thigh_
        << tab << jt.Tcommon_;

    os << nl << "    ";

    for
    (
        register label coefLabel=0;
        coefLabel<janafThermo<equationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        os << jt.highCpCoeffs_[coefLabel] << ' ';
    }

    os << nl << "    ";

    for
    (
        register label coefLabel=0;
        coefLabel<janafThermo<equationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        os << jt.lowCpCoeffs_[coefLabel] << ' ';
    }

    os << endl;

    os.check
    (
        "operator<<(Ostream& os, const janafThermo<equationOfState>& jt)"
    );

    return os;
}


// ************************************************************************* //
