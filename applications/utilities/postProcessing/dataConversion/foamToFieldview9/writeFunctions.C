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

Description

\*---------------------------------------------------------------------------*/

#include "writeFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// C++ version of fwrite_str80 from FieldView/uns/write_binary_uns.c
// Write padded string of 80 char.
bool writeStr80(std::ostream& os, const string& str)
{
    char cBuf[80];

    memset(cBuf, '\0', 80);

    int len = str.size();

    strncpy(cBuf, str.c_str(), (len < 80 ? len : 80));

    os.write(cBuf, 80*sizeof(char));

    return os.good();
}


// Write single integer
bool writeInt(std::ostream& os, int val1)
{
    os.write(reinterpret_cast<char*>(&val1), sizeof(int));

    return os.good();
}


// Write single float
bool writeFloat(std::ostream& os, scalar val1)
{
    float floatVal = val1;

    os.write(reinterpret_cast<char*>(&floatVal), sizeof(float));

    return os.good();
}


// Debug: write raw bytes
void writeBytes(char* start, int nBytes)
{
    cout.setf(std::ios::hex, std::ios::basefield);

    cout<< start << " : ";

    for(int i = 0; i < nBytes; i++)
    {
        cout<< " " << start[i];
    }
    cout << std::endl;

    cout.setf(std::ios::dec);
}


// Debug: write wall flags data
void writeWallFlags(Ostream& os, label cellI, const labelList& wallFlags)
{
    os << "cell " << cellI << " wallsFlags:";
    forAll(wallFlags, wallFaceI)
    {
        os << wallFlags[wallFaceI] << " ";
    }
    os << endl;
}


// ************************************************************************* //
