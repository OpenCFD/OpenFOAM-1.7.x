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

#include "chemistrySolver.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::autoPtr<Foam::chemistrySolver<CompType, ThermoType> >
Foam::chemistrySolver<CompType, ThermoType>::New
(
    ODEChemistryModel<CompType, ThermoType>& model,
    const word& compTypeName,
    const word& thermoTypeName
)
{
    word modelName(model.lookup("chemistrySolver"));

    word chemistrySolverType =
        modelName + '<' + compTypeName + ',' + thermoTypeName + '>';

    Info<< "Selecting chemistrySolver " << modelName << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chemistrySolverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        wordList models = dictionaryConstructorTablePtr_->sortedToc();
        forAll(models, i)
        {
            models[i] = models[i].replace
            (
                '<' + compTypeName + ',' + thermoTypeName + '>',
                ""
            );
        }

        FatalErrorIn
        (
            "chemistrySolver::New"
            "("
                "const ODEChemistryModel&, "
                "const word&, "
                "const word&"
            ")"
        )   << "Unknown chemistrySolver type " << modelName
            << nl << nl << "Valid chemistrySolver types are:" << nl
            << models << nl << exit(FatalError);
    }

    return autoPtr<chemistrySolver<CompType, ThermoType> >
        (cstrIter()(model, modelName));
}


// ************************************************************************* //
