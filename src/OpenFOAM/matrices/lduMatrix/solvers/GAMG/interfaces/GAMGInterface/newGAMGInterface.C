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

#include "GAMGInterface.H"
#include "lduMatrix.H"


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::GAMGInterface> Foam::GAMGInterface::New
(
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
{
    word coupleType(fineInterface.type());

    lduInterfaceConstructorTable::iterator cstrIter =
        lduInterfaceConstructorTablePtr_->find(coupleType);

    if (cstrIter == lduInterfaceConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "GAMGInterface::New"
            "(const lduInterface& fineInterface, "
            "const labelField& localRestrictAddressing, "
            "const labelField& neighbourRestrictAddressing)"
        )   << "Unknown GAMGInterface type " << coupleType << ".\n"
            << "Valid GAMGInterface types are :"
            << lduInterfaceConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<GAMGInterface>
    (
        cstrIter()
        (
            fineInterface,
            localRestrictAddressing,
            neighbourRestrictAddressing
        )
    );
}


// ************************************************************************* //
