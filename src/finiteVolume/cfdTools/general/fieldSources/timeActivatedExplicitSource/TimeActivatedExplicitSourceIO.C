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

#include "TimeActivatedExplicitSource.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TimeActivatedExplicitSource<Type>::writeData(Ostream& os) const
{
    os  << indent << name_ << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    os.writeKeyword("active") << active_ << token::END_STATEMENT << nl;
    os.writeKeyword("timeStart") << timeStart_ << token::END_STATEMENT << nl;
    os.writeKeyword("duration") << duration_ << token::END_STATEMENT << nl;
    os.writeKeyword("selectionMode") << selectionModeTypeToWord(selectionMode_)
        << token::END_STATEMENT << nl;
    os.writeKeyword("volumeMode") << volumeModeTypeToWord(volumeMode_)
        << token::END_STATEMENT << nl;

    if (fieldIds_.size() == 1)
    {
        os.writeKeyword("fieldData") << fieldData_[0].second()
            << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeKeyword("fieldData") << fieldData_ << nl;
    }

    switch (selectionMode_)
    {
        case smPoints:
        {
            os.writeKeyword("points") << nl << indent << points_
                << token::END_STATEMENT << nl;
            break;
        }
        case smCellSet:
        {
            os.writeKeyword("cellSet") << cellSetName_
                << token::END_STATEMENT << nl;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "TimeActivatedExplicitSource<Type>::writeData"
                "("
                    "Ostream&, "
                    "bool"
                ") const"
            )   << "Unknown selectionMode "
                << selectionModeTypeToWord(selectionMode_)
                << abort(FatalError);
        }
    }

    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TimeActivatedExplicitSource<Type>& source
)
{
    source.writeData(os);
    return os;
}


// ************************************************************************* //
