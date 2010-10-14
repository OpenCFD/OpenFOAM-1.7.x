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

#include "incompressible.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::incompressible::incompressible
(
    const specie& sp,
    const scalar rho
)
:
    specie(sp),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::incompressible::incompressible
(
    const word& name,
    const incompressible& ico
)
:
    specie(name, ico),
    rho_(ico.rho_)
{}

inline Foam::autoPtr<Foam::incompressible>
Foam::incompressible::clone() const
{
    return autoPtr<incompressible>(new incompressible(*this));
}

inline Foam::autoPtr<Foam::incompressible>
Foam::incompressible::New(Istream& is)
{
    return autoPtr<incompressible>(new incompressible(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::incompressible::rho(scalar p, scalar T) const
{
    return rho_;
}

inline Foam::scalar Foam::incompressible::psi(scalar, scalar T) const
{
    return 0.0;
}

inline Foam::scalar Foam::incompressible::Z(scalar, scalar) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::incompressible::operator+=(const incompressible& ico)
{
    scalar molr1 = this->nMoles();

    specie::operator+=(ico);

    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ + molr2*ico.rho_;
}

inline void Foam::incompressible::operator-=(const incompressible& ico)
{
    scalar molr1 = this->nMoles();

    specie::operator-=(ico);

    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ - molr2*ico.rho_;
}

inline void Foam::incompressible::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::incompressible Foam::operator+
(
    const incompressible& ico1,
    const incompressible& ico2
)
{
    scalar nMoles = ico1.nMoles() + ico2.nMoles();
    scalar molr1 = ico1.nMoles()/nMoles;
    scalar molr2 = ico2.nMoles()/nMoles;

    return incompressible
    (
        static_cast<const specie&>(ico1)
      + static_cast<const specie&>(ico2),
        molr1*ico1.rho_ + molr2*ico2.rho_
    );
}

inline Foam::incompressible Foam::operator-
(
    const incompressible& ico1,
    const incompressible& ico2
)
{
    scalar nMoles = ico1.nMoles() + ico2.nMoles();
    scalar molr1 = ico1.nMoles()/nMoles;
    scalar molr2 = ico2.nMoles()/nMoles;

    return incompressible
    (
        static_cast<const specie&>(ico1)
      - static_cast<const specie&>(ico2),
        molr1*ico1.rho_ - molr2*ico2.rho_
    );
}

inline Foam::incompressible Foam::operator*
(
    const scalar s,
    const incompressible& ico
)
{
    return incompressible(s*static_cast<const specie&>(ico), ico.rho_);
}

inline Foam::incompressible Foam::operator==
(
    const incompressible& ico1,
    const incompressible& ico2
)
{
    return ico2 - ico1;
}


// ************************************************************************* //
