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

\*---------------------------------------------------------------------------*/

#include "ReactingParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class ParcelType>
Foam::string Foam::ReactingParcel<ParcelType>::propHeader =
    ThermoParcel<ParcelType>::propHeader
  + " mass0"
  + " nPhases(Y1..YN)";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    ThermoParcel<ParcelType>(cloud, is, readFields),
    mass0_(0.0),
    Y_(0),
    pc_(0.0)
{
    if (readFields)
    {
        const ReactingCloud<ParcelType>& cR =
            dynamic_cast<const ReactingCloud<ParcelType>&>(cloud);

        const label nMixture = cR.composition().phaseTypes().size();
        Y_.setSize(nMixture);

        if (is.format() == IOstream::ASCII)
        {
            is >> mass0_ >> Y_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&mass0_),
              + sizeof(mass0_)
            );
            is >> Y_;
        }
    }

    // Check state of Istream
    is.check
    (
        "ReactingParcel<ParcelType>::ReactingParcel"
        "("
            "const Cloud<ParcelType>&, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
void Foam::ReactingParcel<ParcelType>::readFields(Cloud<ParcelType>& cIn)
{
    if (!cIn.size())
    {
        return;
    }

    ReactingCloud<ParcelType>& c =
        dynamic_cast<ReactingCloud<ParcelType>&>(cIn);

    ThermoParcel<ParcelType>::readFields(c);

    IOField<scalar> mass0(c.fieldIOobject("mass0", IOobject::MUST_READ));
    c.checkFieldIOobject(c, mass0);

    label i = 0;
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ReactingParcel<ParcelType>& p = iter();
        p.mass0_ = mass0[i++];
    }

    // Get names and sizes for each Y...
    const wordList& phaseTypes = c.composition().phaseTypes();
    const label nPhases = phaseTypes.size();
    wordList stateLabels(nPhases, "");
    if (c.composition().nPhase() == 1)
    {
        stateLabels = c.composition().stateLabels();
    }


    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ReactingParcel<ParcelType>& p = iter();
        p.Y_.setSize(nPhases, 0.0);
    }

    // Populate Y for each parcel
    forAll(phaseTypes, j)
    {
        IOField<scalar> Y
        (
            c.fieldIOobject
            (
                "Y" + phaseTypes[j] + stateLabels[j],
                 IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter(typename Cloud<ParcelType>, c, iter)
        {
            ReactingParcel<ParcelType>& p = iter();
            p.Y_[j] = Y[i++];
        }
    }
}


template<class ParcelType>
void Foam::ReactingParcel<ParcelType>::writeFields
(
    const Cloud<ParcelType>& cIn
)
{
    const ReactingCloud<ParcelType>& c =
        dynamic_cast<const ReactingCloud<ParcelType>&>(cIn);

    ThermoParcel<ParcelType>::writeFields(c);

    const label np = c.size();

    if (np > 0)
    {
        IOField<scalar> mass0(c.fieldIOobject("mass0", IOobject::NO_READ), np);

        label i = 0;
        forAllConstIter(typename Cloud<ParcelType>, c, iter)
        {
            const ReactingParcel<ParcelType>& p = iter();
            mass0[i++] = p.mass0_;
        }
        mass0.write();

        // Write the composition fractions
        const wordList& phaseTypes = c.composition().phaseTypes();
        wordList stateLabels(phaseTypes.size(), "");
        if (c.composition().nPhase() == 1)
        {
            stateLabels = c.composition().stateLabels();
        }

        forAll(phaseTypes, j)
        {
            IOField<scalar> Y
            (
                c.fieldIOobject
                (
                    "Y" + phaseTypes[j] + stateLabels[j],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter(typename Cloud<ParcelType>, c, iter)
            {
                const ReactingParcel<ParcelType>& p0 = iter();
                Y[i++] = p0.Y()[j];
            }

            Y.write();
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ThermoParcel<ParcelType>&>(p)
            << token::SPACE << p.mass0()
            << token::SPACE << p.Y();
    }
    else
    {
        os  << static_cast<const ThermoParcel<ParcelType>&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.mass0_),
            sizeof(p.mass0())
        );
        os  << p.Y();
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const ReactingParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
