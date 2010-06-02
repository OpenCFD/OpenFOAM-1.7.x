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

Class
    Foam::substance

Description

SourceFiles
    substanceI.H
    substance.C
    substanceIO.C

\*---------------------------------------------------------------------------*/

#ifndef substance_H
#define substance_H

#include "scalar.H"
#include "word.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class substance Declaration
\*---------------------------------------------------------------------------*/

class substance
{
    // Private data

        word name_;
        scalar volFrac_;


public:

    // Constructors

        //- Construct null
        substance()
        {}


    // Member Functions

        // Access

            const word& name() const
            {
                return name_;
            }

            scalar volFrac() const
            {
                return volFrac_;
            }

            bool operator==(const substance& s) const
            {
                return name_ == s.name_ && volFrac_ == s.volFrac_;
            }

            bool operator!=(const substance& s) const
            {
                return !operator==(s);
            }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, substance& s)
        {
            is >> s.name_ >> s.volFrac_;
            return is;
        }

        friend Ostream& operator<<(Ostream& os, const substance& s)
        {
            os << s.name_ << token::SPACE << s.volFrac_;
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
