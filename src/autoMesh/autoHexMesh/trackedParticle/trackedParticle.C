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

\*----------------------------------------------------------------------------*/

#include "trackedParticle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
Foam::trackedParticle::trackedParticle
(
    const Cloud<trackedParticle>& c,
    const vector& position,
    const label celli,
    const point& end,
    const label level,
    const label i,
    const label j
)
:
    ExactParticle<trackedParticle>(c, position, celli),
    end_(end),
    level_(level),
    i_(i),
    j_(j)
{}


//- Construct from Istream
Foam::trackedParticle::trackedParticle
(
    const Cloud<trackedParticle>& c,
    Istream& is,
    bool readFields
)
:
    ExactParticle<trackedParticle>(c, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> end_;
            level_ = readLabel(is);
            i_ = readLabel(is);
            j_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&end_),
                sizeof(end_) + sizeof(level_) + sizeof(i_) + sizeof(j_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "trackedParticle::trackedParticle"
        "(const Cloud<trackedParticle>&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::trackedParticle::move(trackedParticle::trackData& td)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    scalar deltaT = cloud().pMesh().time().deltaT().value();
    scalar tEnd = (1.0 - stepFraction())*deltaT;
    scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // mark visited cell with max level.
        td.maxLevel()[cell()] = max(td.maxLevel()[cell()], level_);

        dt *= trackToFace(end_, td);

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/deltaT;
    }

    return td.keepParticle;
}


bool Foam::trackedParticle::hitPatch
(
    const polyPatch&,
    trackedParticle::trackData& td,
    const label patchI
)
{
    return false;
}


bool Foam::trackedParticle::hitPatch
(
    const polyPatch&,
    int&,
    const label
)
{
    return false;
}


void Foam::trackedParticle::hitWedgePatch
(
    const wedgePolyPatch&,
    trackedParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitWedgePatch
(
    const wedgePolyPatch&,
    int&
)
{}


void Foam::trackedParticle::hitSymmetryPatch
(
    const symmetryPolyPatch&,
    trackedParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitSymmetryPatch
(
    const symmetryPolyPatch&,
    int&
)
{}


void Foam::trackedParticle::hitCyclicPatch
(
    const cyclicPolyPatch&,
    trackedParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitCyclicPatch
(
    const cyclicPolyPatch&,
    int&
)
{}


void Foam::trackedParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackedParticle::trackData& td
)
{
    // Remove particle
    td.switchProcessor = true;
}


void Foam::trackedParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


void Foam::trackedParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackedParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    int&
)
{}


void Foam::trackedParticle::hitPatch
(
    const polyPatch& wpp,
    trackedParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitPatch
(
    const polyPatch& wpp,
    int&
)
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const trackedParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ExactParticle<trackedParticle>&>(p)
            << token::SPACE << p.end_
            << token::SPACE << p.level_
            << token::SPACE << p.i_
            << token::SPACE << p.j_;
    }
    else
    {
        os  << static_cast<const ExactParticle<trackedParticle>&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.end_),
            sizeof(p.end_) + sizeof(p.level_) + sizeof(p.i_) + sizeof(p.j_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const trackedParticle&)");

    return os;
}


// ************************************************************************* //
