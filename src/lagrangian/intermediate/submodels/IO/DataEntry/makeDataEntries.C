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

#include "DataEntry.H"
#include "Constant.H"
#include "Table.H"

#include "label.H"
#include "scalar.H"
#include "vector.H"

namespace Foam
{
    makeDataEntry(label);
    makeDataEntryType(Constant, label);
    makeDataEntryType(Table, label);

    makeDataEntry(scalar);
    makeDataEntryType(Constant, scalar);
    makeDataEntryType(Table, scalar);

    makeDataEntry(vector);
    makeDataEntryType(Constant, vector);
    makeDataEntryType(Table, vector);
};


// ************************************************************************* //
