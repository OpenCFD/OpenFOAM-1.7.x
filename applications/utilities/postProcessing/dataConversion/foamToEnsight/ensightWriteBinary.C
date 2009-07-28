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

#include "ensightWriteBinary.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeEnsDataBinary
(
    const char* val,
    std::ofstream& ensFile
)
{
    char buffer[80] = {0};
    strcpy(buffer, val);
    ensFile.write(buffer, 80*sizeof(char));
}


void writeEnsDataBinary
(
    const int val,
    std::ofstream& ensFile
)
{
    ensFile.write(reinterpret_cast<const char*>(&val), sizeof(int));
}


void writeEnsDataBinary
(
    const scalarField& sf,
    std::ofstream& ensightFile
)
{
    if (sf.size())
    {
        List<float> temp(sf.size());

        forAll(sf, i)
        {
            temp[i] = float(sf[i]);
        }

        ensightFile.write
        (
            reinterpret_cast<char*>(temp.begin()),
            sf.size()*sizeof(float)
        );
    }
}


// ************************************************************************* //

