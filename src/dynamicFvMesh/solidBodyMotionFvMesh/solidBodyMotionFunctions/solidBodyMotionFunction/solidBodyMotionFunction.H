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

Namespace
    Foam::solidBodyMotionFunctions

Description
    Namespace for solid-body motions


Class
    Foam::solidBodyMotionFunction

Description
    Base class for defining solid-body motions

SourceFiles
    solidBodyMotionFunction.C
    newDynamicFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyMotionFunction_H
#define solidBodyMotionFunction_H

#include "Time.H"
#include "dictionary.H"
#include "septernion.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidBodyMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class solidBodyMotionFunction
{
protected:

    // Protected data

        dictionary SBMFCoeffs_;
        const Time& time_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidBodyMotionFunction(const solidBodyMotionFunction&);

        //- Disallow default bitwise assignment
        void operator=(const solidBodyMotionFunction&);


public:

    //- Runtime type information
    TypeName("solidBodyMotionFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidBodyMotionFunction,
            dictionary,
            (const dictionary& SBMFCoeffs, const Time& runTime),
            (SBMFCoeffs, runTime)
        );


    // Constructors

        //- Construct from the SBMFCoeffs dictionary and Time
        solidBodyMotionFunction
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    // Selectors

        //- Select constructed from the SBMFCoeffs dictionary and Time
        static autoPtr<solidBodyMotionFunction> New
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    // Destructor

        virtual ~solidBodyMotionFunction();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
