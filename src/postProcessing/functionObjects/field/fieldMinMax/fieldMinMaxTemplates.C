/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "fieldMinMax.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::fieldMinMax::calcMinMaxFields(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        if (Pstream::master())
        {
            const fieldType& field = obr_.lookupObject<fieldType>(fieldName);
            switch (mode_)
            {
                case mdMag:
                {
                    scalar minValue = min(mag(field)).value();
                    scalar maxValue = max(mag(field)).value();

                    fieldMinMaxFilePtr_() << obr_.time().value() << tab
                        << fieldName << tab << minValue << tab << maxValue
                        << endl;

                    if (log_)
                    {
                        Info<< "fieldMinMax output:" << nl
                            << "    min(mag(" << fieldName << ")) = "
                            << minValue << nl
                            << "    max(mag(" << fieldName << ")) = "
                            << maxValue << nl
                            << endl;
                    }
                    break;
                }
                case mdCmpt:
                {
                    Type minValue = min(field).value();
                    Type maxValue = max(field).value();

                    fieldMinMaxFilePtr_() << obr_.time().value() << tab
                        << fieldName << tab << minValue << tab << maxValue
                        << endl;

                    if (log_)
                    {
                        Info<< "fieldMinMax output:" << nl
                            << "    cmptMin(" << fieldName << ") = "
                            << minValue << nl
                            << "    cmptMax(" << fieldName << ") = "
                            << maxValue << nl
                            << endl;
                    }
                    break;
                }
                default:
                {
                    FatalErrorIn
                    (
                        "Foam::fieldMinMax::calcMinMaxFields"
                        "(const word& fieldName)"
                    )<< "Unknown min/max mode: " << modeTypeNames_[mode_]
                     << exit(FatalError);
                }
            }
        }
    }
}


// ************************************************************************* //
