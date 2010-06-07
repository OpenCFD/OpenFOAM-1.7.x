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

#include "DragModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DragModel<CloudType>::DragModel
(
    const dictionary& dict,
    CloudType& owner
)
:
    dict_(dict),
    owner_(owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DragModel<CloudType>::~DragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::DragModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::DragModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
Foam::scalar Foam::DragModel<CloudType>::utc
(
    const scalar Re,
    const scalar d,
    const scalar mu
) const
{
    const scalar Cd = this->Cd(Re);

    return Cd*Re/d*mu/8.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewDragModel.C"

// ************************************************************************* //

