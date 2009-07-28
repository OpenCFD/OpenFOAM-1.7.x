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
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<radiationModel> radiationModel::New
(
    const volScalarField& T
)
{
    word radiationModelTypeName;

    // Enclose the creation of the radiationPropertiesDict to ensure it is
    // deleted before the radiationModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary radiationPropertiesDict
        (
            IOobject
            (
                "radiationProperties",
                T.mesh().time().constant(),
                T.mesh().objectRegistry::db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        radiationPropertiesDict.lookup("radiationModel")
            >> radiationModelTypeName;
    }

    Info<< "Selecting radiationModel " << radiationModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(radiationModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "radiationModel::New(const volScalarField&)"
        )   << "Unknown radiationModel type " << radiationModelTypeName
            << nl << nl
            << "Valid radiationModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<radiationModel>(cstrIter()(T));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End radiation
} // End namespace Foam

// ************************************************************************* //
