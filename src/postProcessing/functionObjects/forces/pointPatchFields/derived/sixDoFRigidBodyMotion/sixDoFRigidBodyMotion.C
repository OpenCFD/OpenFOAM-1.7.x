/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    FITNESS FOR A PARTICLUAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::applyRestraints()
{
    if (restraints_.empty())
    {
        return;
    }

    if (Pstream::master())
    {
        forAll(restraints_, rI)
        {
            if (report_)
            {
                Info<< "Restraint " << restraintNames_[rI] << ": ";
            }

            // restraint position
            point rP = vector::zero;

            // restraint force
            vector rF = vector::zero;

            // restraint moment
            vector rM = vector::zero;

            restraints_[rI].restrain(*this, rP, rF, rM);

            a() += rF/mass_;

            // Moments are returned in global axes, transforming to
            // body local to add to torque.
            tau() += Q().T() & (rM + ((rP - centreOfMass()) ^ rF));
        }
    }
}


void Foam::sixDoFRigidBodyMotion::applyConstraints(scalar deltaT)
{
    if (constraints_.empty())
    {
        return;
    }

    if (Pstream::master())
    {
        label iteration = 0;

        bool allConverged = true;

        // constraint force accumulator
        vector cFA = vector::zero;

        // constraint moment accumulator
        vector cMA = vector::zero;

        do
        {
            allConverged = true;

            forAll(constraints_, cI)
            {
                if (sixDoFRigidBodyMotionConstraint::debug)
                {
                    Info<< "Constraint " << constraintNames_[cI] << ": ";
                }

                // constraint position
                point cP = vector::zero;

                // constraint force
                vector cF = vector::zero;

                // constraint moment
                vector cM = vector::zero;

                bool constraintConverged = constraints_[cI].constrain
                (
                    *this,
                    cFA,
                    cMA,
                    deltaT,
                    cP,
                    cF,
                    cM
                );

                allConverged = allConverged && constraintConverged;

                // Accumulate constraint force
                cFA += cF;

                // Accumulate constraint moment
                cMA += cM + ((cP - centreOfMass()) ^ cF);
            }

        } while(++iteration < maxConstraintIterations_ && !allConverged);

        if (iteration >= maxConstraintIterations_)
        {
            FatalErrorIn
            (
                "Foam::sixDoFRigidBodyMotion::applyConstraints"
                "(scalar deltaT)"
            )
                << nl << "Maximum number of sixDoFRigidBodyMotion constraint "
                << "iterations ("
                << maxConstraintIterations_
                << ") exceeded." << nl
                << exit(FatalError);
        }

        Info<< "sixDoFRigidBodyMotion constraints converged in "
            << iteration << " iterations" << endl;

        if (report_)
        {
            Info<< "Constraint force: " << cFA << nl
                << "Constraint moment: " << cMA
                << endl;
        }

        // Add the constraint forces and moments to the motion state variables
        a() += cFA/mass_;

        // The moment of constraint forces has already been added
        // during accumulation.  Moments are returned in global axes,
        // transforming to body local
        tau() += Q().T() & cMA;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion()
:
    motionState_(),
    restraints_(),
    restraintNames_(),
    constraints_(),
    constraintNames_(),
    maxConstraintIterations_(0),
    initialCentreOfMass_(vector::zero),
    initialQ_(I),
    momentOfInertia_(diagTensor::one*VSMALL),
    mass_(VSMALL),
    report_(false)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const point& centreOfMass,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau,
    scalar mass,
    const point& initialCentreOfMass,
    const tensor& initialQ,
    const diagTensor& momentOfInertia,
    bool report
)
:
    motionState_
    (
        centreOfMass,
        Q,
        v,
        a,
        pi,
        tau
    ),
    restraints_(),
    restraintNames_(),
    constraints_(),
    constraintNames_(),
    maxConstraintIterations_(0),
    initialCentreOfMass_(initialCentreOfMass),
    initialQ_(initialQ),
    momentOfInertia_(momentOfInertia),
    mass_(mass),
    report_(report)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion(const dictionary& dict)
:
    motionState_(dict),
    restraints_(),
    restraintNames_(),
    constraints_(),
    constraintNames_(),
    maxConstraintIterations_(0),
    initialCentreOfMass_
    (
        dict.lookupOrDefault("initialCentreOfMass", centreOfMass())
    ),
    initialQ_
    (
        dict.lookupOrDefault("initialOrientation", Q())
    ),
    momentOfInertia_(dict.lookup("momentOfInertia")),
    mass_(readScalar(dict.lookup("mass"))),
    report_(dict.lookupOrDefault<Switch>("report", false))
{
    addRestraints(dict);

    addConstraints(dict);
}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const sixDoFRigidBodyMotion& sDoFRBM
)
:
    motionState_(sDoFRBM.motionState()),
    restraints_(sDoFRBM.restraints()),
    restraintNames_(sDoFRBM.restraintNames()),
    constraints_(sDoFRBM.constraints()),
    constraintNames_(sDoFRBM.constraintNames()),
    maxConstraintIterations_(sDoFRBM.maxConstraintIterations()),
    initialCentreOfMass_(sDoFRBM.initialCentreOfMass()),
    initialQ_(sDoFRBM.initialQ()),
    momentOfInertia_(sDoFRBM.momentOfInertia()),
    mass_(sDoFRBM.mass()),
    report_(sDoFRBM.report())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::~sixDoFRigidBodyMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::addRestraints
(
    const dictionary& dict
)
{
    if (dict.found("restraints"))
    {
        const dictionary& restraintDict = dict.subDict("restraints");

        label i = 0;

        restraints_.setSize(restraintDict.size());

        restraintNames_.setSize(restraintDict.size());

        forAllConstIter(IDLList<entry>, restraintDict, iter)
        {
            if (iter().isDict())
            {
                // Info<< "Adding restraint: " << iter().keyword() << endl;

                restraints_.set
                (
                    i,
                    sixDoFRigidBodyMotionRestraint::New(iter().dict())
                );

                restraintNames_[i] = iter().keyword();

                i++;
            }
        }

        restraints_.setSize(i);

        restraintNames_.setSize(i);
    }
}


void Foam::sixDoFRigidBodyMotion::addConstraints
(
    const dictionary& dict
)
{
    if (dict.found("constraints"))
    {
        const dictionary& constraintDict = dict.subDict("constraints");

        label i = 0;

        constraints_.setSize(constraintDict.size());

        constraintNames_.setSize(constraintDict.size());

        forAllConstIter(IDLList<entry>, constraintDict, iter)
        {
            if (iter().isDict())
            {
                // Info<< "Adding constraint: " << iter().keyword() << endl;

                constraints_.set
                (
                    i,
                    sixDoFRigidBodyMotionConstraint::New(iter().dict())
                );

                constraintNames_[i] = iter().keyword();

                i++;
            }
        }

        constraints_.setSize(i);

        constraintNames_.setSize(i);

        if (!constraints_.empty())
        {
            maxConstraintIterations_ = readLabel
            (
                constraintDict.lookup("maxIterations")
            );
        }
    }
}


void Foam::sixDoFRigidBodyMotion::updatePosition
(
    scalar deltaT
)
{
    // First leapfrog velocity adjust and motion part, required before
    // force calculation

    if (Pstream::master())
    {
        v() += 0.5*deltaT*a();

        pi() += 0.5*deltaT*tau();

        // Leapfrog move part
        centreOfMass() += deltaT*v();

        // Leapfrog orientation adjustment

        rotate(Q(), pi(), deltaT);
    }

    Pstream::scatter(motionState_);
}


void Foam::sixDoFRigidBodyMotion::updateForce
(
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT
)
{
    // Second leapfrog velocity adjust part, required after motion and
    // force calculation

    if (Pstream::master())
    {
        a() = fGlobal/mass_;

        tau() = (Q().T() & tauGlobal);

        applyRestraints();

        applyConstraints(deltaT);

        v() += 0.5*deltaT*a();

        pi() += 0.5*deltaT*tau();

        if(report_)
        {
            status();
        }
    }

    Pstream::scatter(motionState_);
}


void Foam::sixDoFRigidBodyMotion::updateForce
(
    const pointField& positions,
    const vectorField& forces,
    scalar deltaT
)
{
    vector a = vector::zero;

    vector tau = vector::zero;

    if (Pstream::master())
    {
        forAll(positions, i)
        {
            const vector& f = forces[i];

            a += f/mass_;

            tau += Q().T() & ((positions[i] - centreOfMass()) ^ f);
        }
    }

    updateForce(a, tau, deltaT);
}


Foam::point Foam::sixDoFRigidBodyMotion::predictedPosition
(
    const point& pInitial,
    const vector& deltaForce,
    const vector& deltaMoment,
    scalar deltaT
) const
{
    vector vTemp = v() + deltaT*(a() + deltaForce/mass_);

    vector piTemp = pi() + deltaT*(tau() + (Q().T() & deltaMoment));

    point centreOfMassTemp = centreOfMass() + deltaT*vTemp;

    tensor QTemp = Q();

    rotate(QTemp, piTemp, deltaT);

    return
    (
        centreOfMassTemp
      + (QTemp & initialQ_.T() & (pInitial - initialCentreOfMass_))
    );
}


Foam::vector Foam::sixDoFRigidBodyMotion::predictedOrientation
(
    const vector& vInitial,
    const vector& deltaMoment,
    scalar deltaT
) const
{
    vector piTemp = pi() + deltaT*(tau() + (Q().T() & deltaMoment));

    tensor QTemp = Q();

    rotate(QTemp, piTemp, deltaT);

    return (QTemp & initialQ_.T() & vInitial);
}


void Foam::sixDoFRigidBodyMotion::status() const
{
    Info<< "Centre of mass: " << centreOfMass() << nl
        << "Linear velocity: " << v() << nl
        << "Angular velocity: " << omega()
        << endl;
}


// ************************************************************************* //
