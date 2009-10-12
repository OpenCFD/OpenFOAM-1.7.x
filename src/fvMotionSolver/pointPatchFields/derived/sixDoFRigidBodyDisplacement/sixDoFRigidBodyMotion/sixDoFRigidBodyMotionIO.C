/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "sixDoFRigidBodyMotion.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::write(Ostream& os) const
{
    motionState_.write(os);

    os.writeKeyword("refCentreOfMass")
        << refCentreOfMass_ << token::END_STATEMENT << nl;
    os.writeKeyword("momentOfInertia")
        << momentOfInertia_ << token::END_STATEMENT << nl;
    os.writeKeyword("mass")
        << mass_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, sixDoFRigidBodyMotion& sDoFRBM)
{
    is  >> sDoFRBM.motionState_
        >> sDoFRBM.refCentreOfMass_
        >> sDoFRBM.momentOfInertia_
        >> sDoFRBM.mass_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::sixDoFRigidBodyMotion&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sixDoFRigidBodyMotion& sDoFRBM
)
{
    os  << sDoFRBM.motionState()
        << token::SPACE << sDoFRBM.refCentreOfMass()
        << token::SPACE << sDoFRBM.momentOfInertia()
        << token::SPACE << sDoFRBM.mass() ;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::sixDoFRigidBodyMotion&)"
    );

    return os;
}


// ************************************************************************* //
