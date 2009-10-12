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

Class
    Foam::sixDoFRigidBodyMotion

Description
    Six degree of freedom motion for a rigid body.  Angular momentum stored in
    body fixed reference frame.  Reference orientation of the body must align
    with the cartesian axes such that the Inertia tensor is in principle
    component form.

    Symplectic motion as per:

    title = {Symplectic splitting methods for rigid body molecular dynamics},
    publisher = {AIP},
    year = {1997},
    journal = {The Journal of Chemical Physics},
    volume = {107},
    number = {15},
    pages = {5840-5851},
    url = {http://link.aip.org/link/?JCP/107/5840/1},
    doi = {10.1063/1.474310}

SourceFiles
    sixDoFRigidBodyMotionI.H
    sixDoFRigidBodyMotion.C
    sixDoFRigidBodyMotionIO.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotion_H
#define sixDoFRigidBodyMotion_H

#include "sixDoFRigidBodyMotionState.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class sixDoFRigidBodyMotion;
Istream& operator>>(Istream&, sixDoFRigidBodyMotion&);
Ostream& operator<<(Ostream&, const sixDoFRigidBodyMotion&);


/*---------------------------------------------------------------------------*\
                      Class sixDoFRigidBodyMotion Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotion
{
    // Private data

        // state data object
        sixDoFRigidBodyMotionState motionState_;

        //- Centre of mass of reference state
        point refCentreOfMass_;

        //- Moment of inertia of the body in reference configuration
        diagTensor momentOfInertia_;

        //- Mass of the body
        scalar mass_;


    // Private Member Functions

        //- Calculate the rotation tensor around the body reference
        //  frame x-axis by the given angle
        inline tensor rotationTensorX(scalar deltaT) const;

        //- Calculate the rotation tensor around the body reference
        //  frame y-axis by the given angle
        inline tensor rotationTensorY(scalar deltaT) const;

        //- Calculate the rotation tensor around the body reference
        //  frame z-axis by the given angle
        inline tensor rotationTensorZ(scalar deltaT) const;


public:

    // Constructors

        //- Construct null
        sixDoFRigidBodyMotion();

        //- Construct from components
        sixDoFRigidBodyMotion
        (
            const point& centreOfMass,
            const tensor& Q,
            const vector& v,
            const vector& a,
            const vector& pi,
            const vector& tau,
            scalar mass,
            const point& refCentreOfMass,
            const diagTensor& momentOfInertia
        );

        //- Construct from dictionary
        sixDoFRigidBodyMotion(const dictionary& dict);

        //- Construct as copy
        sixDoFRigidBodyMotion(const sixDoFRigidBodyMotion&);


    //- Destructor
    ~sixDoFRigidBodyMotion();


    // Member Functions

        void updatePosition
        (
            scalar deltaT
        );

        void updateForce
        (
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT
        );

        void updateForce
        (
            const pointField& positions,
            const vectorField& forces,
            scalar deltaT
        );

        tmp<pointField> generatePositions(const pointField& pts) const;

        // Access

            //- Return access to the motion state
            inline const sixDoFRigidBodyMotionState& motionState() const;

            //- Return access to the centre of mass
            inline const point& centreOfMass() const;

            //- Return access to the centre of mass
            inline const point& refCentreOfMass() const;

            //- Return access to the inertia tensor
            inline const diagTensor& momentOfInertia() const;

            //- Return access to the mass
            inline scalar mass() const;

            //- Return access to the orientation
            inline const tensor& Q() const;

            //- Return access to velocity
            inline const vector& v() const;

            //- Return access to acceleration
            inline const vector& a() const;

            //- Return access to angular momentum
            inline const vector& pi() const;

            //- Return access to torque
            inline const vector& tau() const;


        // Edit

            //- Return non-const access to the centre of mass
            inline point& centreOfMass();

            //- Return access to the centre of mass
            inline point& refCentreOfMass();

            //- Return non-const access to the inertia tensor
            inline diagTensor& momentOfInertia();

            //- Return non-const access to the mass
            inline scalar& mass();

            //- Return non-const access to the orientation
            inline tensor& Q();

            //- Return non-const access to vector
            inline vector& v();

            //- Return non-const access to acceleration
            inline vector& a();

            //- Return non-const access to angular momentum
            inline vector& pi();

            //- Return non-const access to torque
            inline vector& tau();


        //- Write
        void write(Ostream&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, sixDoFRigidBodyMotion&);
        friend Ostream& operator<<(Ostream&, const sixDoFRigidBodyMotion&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDoFRigidBodyMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
