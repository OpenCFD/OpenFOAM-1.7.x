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

#include "injectorType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(injectorType, 0);
defineRunTimeSelectionTable(injectorType, dictionary);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::injectorType::injectorType
(
    const Foam::Time&,
    const Foam::dictionary&
)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::injectorType> Foam::injectorType::New
(
    const Time& t,
    const dictionary& dict
)
{
    word injectorTypeName
    (
        dict.lookup("injectorType")
    );

    Info<< "Selecting injectorType "
         << injectorTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(injectorTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "injectorType::New(const dictionary&) : " << endl
            << "    unknown injectorType type "
            << injectorTypeName
            << ", constructor not in hash table" << endl << endl
            << "    Valid injector types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<injectorType>(cstrIter()(t, dict));

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::injectorType::~injectorType()
{}

Foam::scalar Foam::injectorType::getTableValue
(
    const List<pair>& table,
    const scalar value
) const
{
    // iterator
    label i = 0;

    // max items
    label maxRow = table.size() - 1;

    // check lower bound
    if (value < table[0][0])
    {
        return table[0][1];
    }
    // check upper bound
    else if (value > table[maxRow][0])
    {
        return table[maxRow][1];
    }
    // interpolate intermediate value
    else
    {
        while
        (
            (i < maxRow-1) && (table[i+1][0] < value)
        )
        {
            i++;
        }
        // value sits bewteen table[i][0] and table[i+1][0]
        return table[i][1] 
               + (value-table[i][0])/(table[i+1][0]-table[i][0])
               * (table[i+1][1]-table[i][1]);
    }
}

Foam::scalar Foam::injectorType::integrateTable
(
    const List<pair>& table,
    const scalar value
) const
{
    label N = table.size() - 1;
    scalar sum = 0.0;
    scalar t = max(table[0][0], min(value, table[N][0]));

    label i = 0;
    while
    (
        (i < N - 1)
     && (table[i+1][0] < t)
    )
    {
        scalar deltaH = table[i+1][1] + table[i][1];
        scalar deltaT = table[i+1][0] - table[i][0];
        sum += 0.5*deltaH*deltaT;
        i++;
    }

    scalar interpolatedValue =
        table[i][1]
      + (t - table[i][0])
      * (table[i+1][1] - table[i][1])
      / (table[i+1][0] - table[i][0]);

    sum +=
        0.5*(interpolatedValue + table[i][1])
       *(t - table[i][0]);

    return sum;
}

Foam::scalar Foam::injectorType::integrateTable
(
    const List<pair>& table
) const
{
    scalar integratedTable = 0.0;
    for (label i=0; i < table.size() - 1; i++)
    {
        scalar deltaH = table[i+1][1] + table[i][1];
        scalar deltaT = table[i+1][0] - table[i][0];
        integratedTable += 0.5*deltaH*deltaT;
    }

    return integratedTable;
}

// ************************************************************************* //
