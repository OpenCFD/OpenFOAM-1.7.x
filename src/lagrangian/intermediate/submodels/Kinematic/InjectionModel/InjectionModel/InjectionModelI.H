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

#include "InjectionModel.H"

template<class CloudType>
const Foam::dictionary& Foam::InjectionModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const CloudType& Foam::InjectionModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType& Foam::InjectionModel<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::InjectionModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::timeStart() const
{
    return SOI_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::volumeTotal() const
{
    return volumeTotal_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::massTotal() const
{
    return massTotal_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::massInjected() const
{
    return massInjected_;
}


template<class CloudType>
Foam::label Foam::InjectionModel<CloudType>::nInjections() const
{
    return nInjections_;
}


template<class CloudType>
Foam::label Foam::InjectionModel<CloudType>::parcelsAddedTotal() const
{
    return parcelsAddedTotal_;
}


// ************************************************************************* //
