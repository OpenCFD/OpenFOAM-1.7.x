/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::interactionLists::mesh() const
{
    return mesh_;
}


const Foam::directInteractionList& Foam::interactionLists::dil() const
{
    return dil_;
}


inline const Foam::referredCellList& Foam::interactionLists::ril() const
{
    return ril_;
}


inline Foam::referredCellList& Foam::interactionLists::ril()
{
    return ril_;
}


inline const Foam::List<Foam::sendingReferralList>&
    Foam::interactionLists::cellSendingReferralLists() const
{
    return cellSendingReferralLists_;
}


inline const Foam::List<Foam::receivingReferralList>&
    Foam::interactionLists::cellReceivingReferralLists() const
{
    return cellReceivingReferralLists_;
}


inline Foam::label Foam::interactionLists::nInteractingProcs() const
{
    return cellReceivingReferralLists_.size();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
