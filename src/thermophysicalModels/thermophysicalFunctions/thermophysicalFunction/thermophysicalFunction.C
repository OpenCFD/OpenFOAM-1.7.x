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

#include "error.H"

#include "thermophysicalFunction.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermophysicalFunction, 0);
defineRunTimeSelectionTable(thermophysicalFunction, Istream);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<thermophysicalFunction> thermophysicalFunction::New(Istream& is)
{
    if (debug)
    {
        Info<< "thermophysicalFunction::New(Istream&) : "
            << "constructing thermophysicalFunction"
            << endl;
    }

    word thermophysicalFunctionType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(thermophysicalFunctionType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorIn("thermophysicalFunction::New(Istream&)")
            << "Unknown thermophysicalFunction type "
            << thermophysicalFunctionType
            << endl << endl
            << "Valid thermophysicalFunction types are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<thermophysicalFunction>(cstrIter()(is));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
