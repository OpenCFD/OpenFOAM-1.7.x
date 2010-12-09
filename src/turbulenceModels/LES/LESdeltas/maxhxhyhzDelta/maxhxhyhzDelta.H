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

Class
    Foam::maxhxhyhzDelta

Description
    maxhxhyhzDelta takes the maximum of the three dimensions per cell:
    max(hx, hy, hz). Valid for structures hexahedral cells only.


SourceFiles
    maxhxhyhzDelta.C

\*---------------------------------------------------------------------------*/

#ifndef maxhxhyhzDeltaDelta_H
#define maxhxhyhzDeltaDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class maxhxhyhzDelta Declaration
\*---------------------------------------------------------------------------*/

class maxhxhyhzDelta
:
    public LESdelta
{
    // Private data

        scalar deltaCoeff_; //


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        maxhxhyhzDelta(const maxhxhyhzDelta&);
        void operator=(const maxhxhyhzDelta&);

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("maxhxhyhzDelta");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        maxhxhyhzDelta
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        );


    //- Destructor
    virtual ~maxhxhyhzDelta()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
