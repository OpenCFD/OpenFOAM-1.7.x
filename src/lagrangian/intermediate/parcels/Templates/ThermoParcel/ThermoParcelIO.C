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

#include "ThermoParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class ParcelType>
Foam::string Foam::ThermoParcel<ParcelType>::propHeader =
    KinematicParcel<ParcelType>::propHeader
  + " T"
  + " cp";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    KinematicParcel<ParcelType>(cloud, is, readFields),
    T_(0.0),
    cp_(0.0),
    Tc_(0.0),
    cpc_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            T_ = readScalar(is);
            cp_ = readScalar(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&T_),
              + sizeof(T_)
              + sizeof(cp_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "ThermoParcel::ThermoParcel(const Cloud<ParcelType>&, Istream&, bool)"
    );
}


template<class ParcelType>
void Foam::ThermoParcel<ParcelType>::readFields(Cloud<ParcelType>& c)
{
    if (!c.size())
    {
        return;
    }

    KinematicParcel<ParcelType>::readFields(c);

    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ));
    c.checkFieldIOobject(c, T);

    IOField<scalar> cp(c.fieldIOobject("cp", IOobject::MUST_READ));
    c.checkFieldIOobject(c, cp);


    label i = 0;
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ThermoParcel<ParcelType>& p = iter();

        p.T_ = T[i];
        p.cp_ = cp[i];
        i++;
    }
}


template<class ParcelType>
void Foam::ThermoParcel<ParcelType>::writeFields(const Cloud<ParcelType>& c)
{
    KinematicParcel<ParcelType>::writeFields(c);

    label np =  c.size();

    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> cp(c.fieldIOobject("cp", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<ParcelType>, c, iter)
    {
        const ThermoParcel<ParcelType>& p = iter();

        T[i] = p.T_;
        cp[i] = p.cp_;
        i++;
    }

    T.write();
    cp.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ThermoParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const KinematicParcel<ParcelType>&>(p)
            << token::SPACE << p.T()
            << token::SPACE << p.cp();
    }
    else
    {
        os  << static_cast<const KinematicParcel<ParcelType>&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.T_),
            sizeof(p.T()) + sizeof(p.cp())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const ThermoParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
