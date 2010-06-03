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

Class
    Foam::injector

Description
    Spray injector

\*---------------------------------------------------------------------------*/

#ifndef injector_H
#define injector_H

#include "vector.H"
#include "List.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "Time.H"
#include "injectorType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class injector Declaration
\*---------------------------------------------------------------------------*/

class injector
{
    // Private data

        dictionary injectorDict_;

        autoPtr<injectorType> properties_;


public:

    // Constructors
    
        //- Constructor from Istream
        injector(const Time& t, Istream& is);

        //- Return clone
        inline autoPtr<injector> clone() const;

        //- Return a pointer to a new injector created on freestore
        //  from Istream
        class iNew
        {
            const Time& time_;

        public:

            iNew(const Time& t)
            :
                time_(t)
            {}

            autoPtr<injector> operator()(Istream& is) const
            {
                return autoPtr<injector>(new injector(time_, is));
            }
        };


    // Member Functions

        inline autoPtr<injectorType>& properties();
        inline const autoPtr<injectorType>& properties() const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const injector&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "injectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
