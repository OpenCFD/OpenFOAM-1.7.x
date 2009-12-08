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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::tensor
Foam::sixDoFRigidBodyMotion::rotationTensorX(scalar phi) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


inline Foam::tensor
Foam::sixDoFRigidBodyMotion::rotationTensorY(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


inline Foam::tensor
Foam::sixDoFRigidBodyMotion::rotationTensorZ(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


inline const Foam::sixDoFRigidBodyMotionState&
Foam::sixDoFRigidBodyMotion::motionState() const
{
    return motionState_;
}


inline const Foam::point& Foam::sixDoFRigidBodyMotion::centreOfMass() const
{
    return motionState_.centreOfMass();
}


inline const Foam::point& Foam::sixDoFRigidBodyMotion::refCentreOfMass() const
{
    return refCentreOfMass_;
}


inline const Foam::diagTensor&
Foam::sixDoFRigidBodyMotion::momentOfInertia() const
{
    return momentOfInertia_;
}


inline Foam::scalar Foam::sixDoFRigidBodyMotion::mass() const
{
    return mass_;
}


inline const Foam::tensor& Foam::sixDoFRigidBodyMotion::Q() const
{
    return motionState_.Q();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::v() const
{
    return motionState_.v();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::a() const
{
    return motionState_.a();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::pi() const
{
    return motionState_.pi();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::tau() const
{
    return motionState_.tau();
}


inline Foam::point& Foam::sixDoFRigidBodyMotion::centreOfMass()
{
    return motionState_.centreOfMass();
}


inline Foam::point& Foam::sixDoFRigidBodyMotion::refCentreOfMass()
{
    return refCentreOfMass_;
}


inline Foam::diagTensor& Foam::sixDoFRigidBodyMotion::momentOfInertia()
{
    return momentOfInertia_;
}


inline Foam::scalar& Foam::sixDoFRigidBodyMotion::mass()
{
    return mass_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotion::Q()
{
    return motionState_.Q();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotion::v()
{
    return motionState_.v();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotion::a()
{
    return motionState_.a();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotion::pi()
{
    return motionState_.pi();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotion::tau()
{
    return motionState_.tau();
}


// ************************************************************************* //
