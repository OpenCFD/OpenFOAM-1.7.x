/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "Particle.H"
#include "IOstreams.H"
#include "IOPosition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::string Foam::Particle<ParticleType>::propHeader =
    "(Px Py Pz) cellI origProc origId";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Particle<ParticleType>::Particle
(
    const Cloud<ParticleType>& cloud,
    Istream& is,
    bool readFields
)
:
    cloud_(cloud),
    facei_(-1),
    stepFraction_(0.0),
    origProc_(Pstream::myProcNo()),
    origId_(-1)
{

    // readFields : read additional data. Should be consistent with writeFields.

    if (is.format() == IOstream::ASCII)
    {
        is >> position_ >> celli_;
        if (readFields)
        {
            is >> origProc_ >> origId_;
        }
        else
        {
            origId_ = cloud_.getNewParticleID();
        }
    }
    else
    {
        // In binary read all particle data - needed for parallel transfer
        if (readFields)
        {
            is.read
            (
                reinterpret_cast<char*>(&position_),
                sizeof(position_)
              + sizeof(celli_)
              + sizeof(facei_)
              + sizeof(stepFraction_)
              + sizeof(origProc_)
              + sizeof(origId_)
            );
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&position_),
                sizeof(position_)
              + sizeof(celli_)
              + sizeof(facei_)
              + sizeof(stepFraction_)
            );
            origId_ = cloud_.getNewParticleID();
        }
    }

    if (celli_ == -1)
    {
        celli_ = cloud_.pMesh().findCell(position_);
    }

    // Check state of Istream
    is.check("Particle<ParticleType>::Particle(Istream&)");
}


template<class ParticleType>
void Foam::Particle<ParticleType>::readFields
(
    Cloud<ParticleType>& c
)
{
    if (!c.size())
    {
        return;
    }

    IOobject procIO(c.fieldIOobject("origProcId", IOobject::MUST_READ));

    if (procIO.headerOk())
    {
        IOField<label> origProcId(procIO);
        c.checkFieldIOobject(c, origProcId);
        IOField<label> origId(c.fieldIOobject("origId", IOobject::MUST_READ));
        c.checkFieldIOobject(c, origId);

        label i = 0;
        forAllIter(typename Cloud<ParticleType>, c, iter)
        {
            ParticleType& p = iter();

            p.origProc_ = origProcId[i];
            p.origId_ = origId[i];
            i++;
        }
    }
}


template<class ParticleType>
void Foam::Particle<ParticleType>::writeFields
(
    const Cloud<ParticleType>& c
)
{
    // Write the cloud position file
    IOPosition<ParticleType> ioP(c);
    ioP.write();

    label np =  c.size();

    IOField<label> origProc
    (
        c.fieldIOobject
        (
            "origProcId",
            IOobject::NO_READ
        ),
        np
    );
    IOField<label> origId(c.fieldIOobject("origId", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<ParticleType>, c, iter)
    {
        origProc[i] = iter().origProc_;
        origId[i] = iter().origId_;
        i++;
    }

    origProc.write();
    origId.write();
}


template<class ParticleType>
void Foam::Particle<ParticleType>::write(Ostream& os, bool writeFields) const
{
    if (os.format() == IOstream::ASCII)
    {
        if (writeFields)
        {
            // Write the additional entries
            os << position_
               << token::SPACE << celli_
               << token::SPACE << origProc_
               << token::SPACE << origId_;
        }
        else
        {
            os << position_
               << token::SPACE << celli_;
        }
    }
    else
    {
        // In binary write both celli_ and facei_, needed for parallel transfer
        if (writeFields)
        {
            os.write
            (
                reinterpret_cast<const char*>(&position_),
                sizeof(position_)
              + sizeof(celli_)
              + sizeof(facei_)
              + sizeof(stepFraction_)
              + sizeof(origProc_)
              + sizeof(origId_)
            );
        }
        else
        {
            os.write
            (
                reinterpret_cast<const char*>(&position_),
                sizeof(position_)
              + sizeof(celli_)
              + sizeof(facei_)
              + sizeof(stepFraction_)
            );
        }
    }

    // Check state of Ostream
    os.check("Particle<ParticleType>::write(Ostream& os, bool) const");
}


template<class ParticleType>
Foam::Ostream& Foam::operator<<(Ostream& os, const Particle<ParticleType>& p)
{
    // Write all data
    p.write(os, true);

    return os;
}


// ************************************************************************* //
