/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "explicitSource.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::explicitSource::writeData(Ostream& os) const
{
    os  << indent << name_ << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    os.writeKeyword("volumeMode") << volumeModeTypeToWord(volumeMode_)
        << token::END_STATEMENT << nl;

    if (scalarFields_.size() > 0)
    {
        os.writeKeyword("scalarFields") << scalarFields_
            << token::END_STATEMENT << nl;
    }

    if (vectorFields_.size() > 0)
    {
        os.writeKeyword("vectorFields") << vectorFields_
            << token::END_STATEMENT << nl;
    }

    os << decrIndent << indent << token::END_BLOCK << endl;
}


bool Foam::explicitSource::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {
        const dictionary& sourceDict = dict.subDict(name());
        const dictionary& subDictCoeffs = sourceDict.subDict(typeName + "Coeffs");
        setFieldData(subDictCoeffs.subDict("fieldData"));
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const explicitSource& source)
{
    source.writeData(os);
    return os;
}


// ************************************************************************* //
