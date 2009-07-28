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

Class
    Foam::solution

Description
    Selector class for relaxation factors, solver type and solution.

SourceFiles
    solution.C

\*---------------------------------------------------------------------------*/

#ifndef solution_H
#define solution_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class solution Declaration
\*---------------------------------------------------------------------------*/

class solution
:
    public IOdictionary
{
    // Private data

        //- Dictionary of relaxation factors for all the fields
        dictionary relaxationFactors_;

        //- Optional default relaxation factor for all the fields
        scalar defaultRelaxationFactor_;

        //- Dictionary of solver parameters for all the fields
        dictionary solvers_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        solution(const solution&);
        void operator=(const solution&);


public:

    //- Update from older solver controls syntax
    //  Usually verbose, since we want to know about the changes
    //  Returns the number of settings changed
    static label upgradeSolverDict(dictionary& dict, const bool verbose=true);

    //- Debug switch
    static int debug;


    // Constructors

        //- Construct for given objectRegistry and dictionary
        solution(const objectRegistry& obr, const fileName& dictName);


    // Member Functions

        // Access

            //- Return the selected sub-dictionary of solvers if the "select"
            //  keyword is given, otherwise return the complete dictionary
            const dictionary& solutionDict() const;

            //- Return true if the relaxation factor is given for the field
            bool relax(const word& name) const;

            //- Return the relaxation factor for the given field
            scalar relaxationFactor(const word& name) const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solverDict(const word& name) const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solver(const word& name) const;

        // Read

            //- Read the solution dictionary
            bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
