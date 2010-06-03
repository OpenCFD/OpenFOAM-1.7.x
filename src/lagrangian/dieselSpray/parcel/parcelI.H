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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const List<word>& parcel::liquidNames() const
{
    return liquidComponents_;
}

inline const List<word>& parcel::fuelNames() const
{
    return liquidComponents_;
}

inline scalar parcel::d() const
{
    return d_;
}

inline scalar& parcel::d()
{
    return d_;
}

inline scalar parcel::T() const
{
    return T_;
}

inline scalar& parcel::T()
{
    return T_;
}

inline scalar parcel::m() const
{
    return m_;
}

inline scalar& parcel::m()
{
    return m_;
}

inline scalar parcel::dev() const
{
    return y_;
}

inline scalar& parcel::dev()
{
    return y_;
}

inline scalar parcel::ddev() const
{
    return yDot_;
}

inline scalar& parcel::ddev()
{
    return yDot_;
}

inline scalar parcel::ct() const
{
    return ct_;
}

inline scalar& parcel::ct()
{
    return ct_;
}

inline scalar& parcel::ms()
{
    return ms_;
}

inline scalar parcel::ms() const
{
    return ms_;
}

inline scalar& parcel::tTurb()
{
    return tTurb_;
}

inline scalar parcel::tTurb() const
{
    return tTurb_;
}

inline scalar& parcel::liquidCore()
{
    return liquidCore_;
}

inline scalar parcel::liquidCore() const
{
    return liquidCore_;
}

inline scalar& parcel::injector()
{
    return injector_;
}

inline scalar parcel::injector() const
{
    return injector_;
}

inline const vector& parcel::U() const
{
    return U_;
}

inline vector& parcel::U()
{
    return U_;
}

inline const vector& parcel::Uturb() const
{
    return Uturb_;
}

inline vector& parcel::Uturb()
{
    return Uturb_;
}

inline const vector& parcel::n() const
{
    return n_;
}

inline vector& parcel::n()
{
    return n_;
}

inline const scalarField& parcel::X() const
{
    return X_;
}

inline scalarField& parcel::X()
{
    return X_;
}

inline scalar& parcel::tMom()
{
    return tMom_;
}

inline scalar parcel::tMom() const
{
    return tMom_;
}

inline vector parcel::Urel(const vector& U) const
{
    return U_ - U;
}

inline void parcel::correctNormal
(
    const vector& sym
)
{
    scalar d = position() & sym;
    n_ = position() - d*sym;
    n_ /= mag(n_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
