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

#include "Table.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Table<Type>::Table(const word& entryName, Istream& is)
:
    DataEntry<Type>(entryName),
    table_(is)
{
    if (!table_.size())
    {
        FatalErrorIn("Foam::Table<Type>::Table(const Istream&)")
            << "Table for entry " << this->name_ << " is invalid (empty)"
            << nl << exit(FatalError);
    }
}


template<class Type>
Foam::Table<Type>::Table(const Table<Type>& tbl)
:
    DataEntry<Type>(tbl),
    table_(tbl.table_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Table<Type>::~Table()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Table<Type>::value(const scalar x) const
{
    // Return zero if out of bounds
    if (x < table_[0].first() || x > table_[table_.size()-1].first())
    {
        return pTraits<Type>::zero;
    }

    // Find i such that x(i) < x < x(i+1)
    label i = 0;
    while ((table_[i+1].first() < x) && (i+1 < table_.size()))
    {
        i++;
    }

    // Linear interpolation to find value. Note constructor needed for
    // Table<label> to convert intermediate scalar back to label.
    return Type
    (
        (x - table_[i].first())/(table_[i+1].first() - table_[i].first())
      * (table_[i+1].second() - table_[i].second())
      + table_[i].second()
    );
}


template<class Type>
Type Foam::Table<Type>::integrate(const scalar x1, const scalar x2) const
{
    // Initialise return value
    Type sum = pTraits<Type>::zero;

    // Return zero if out of bounds
    if ((x1 > table_[table_.size()-1].first()) || (x2 < table_[0].first()))
    {
        return sum;
    }

    // Find next index greater than x1
    label id1 = 0;
    while ((table_[id1].first() < x1) && (id1 < table_.size()))
    {
        id1++;
    }

    // Find next index less than x2
    label id2 = table_.size() - 1;
    while ((table_[id2].first() > x2) && (id2 >= 1))
    {
        id2--;
    }

    if ((id1 - id2) == 1)
    {
        // x1 and x2 lie within 1 interval
        sum = 0.5*(value(x1) + value(x2))*(x2 - x1);
    }
    else
    {
        // x1 and x2 cross multiple intervals

        // Integrate table body
        for (label i=id1; i<id2; i++)
        {
            sum +=
                (table_[i].second() + table_[i+1].second())
              * (table_[i+1].first() - table_[i].first());
        }
        sum *= 0.5;

        // Add table ends (partial segments)
        if (id1 > 0)
        {
            sum += 0.5
              * (value(x1) + table_[id1].second())
              * (table_[id1].first() - x1);
        }
        if (id2 < table_.size() - 1)
        {
            sum += 0.5
              * (table_[id2].second() + value(x2))
              * (x2 - table_[id2].first());
        }
    }

    return sum;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "TableIO.C"


// ************************************************************************* //
