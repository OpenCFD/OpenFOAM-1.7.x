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
    Foam::basicPsiThermo

Description
    Basic thermodynamic properties based on compressibility

SourceFiles
    basicPsiThermo.C
    newBasicPsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicPsiThermo_H
#define basicPsiThermo_H

#include "basicThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicPsiThermo Declaration
\*---------------------------------------------------------------------------*/

class basicPsiThermo
:
    public basicThermo
{

protected:

    // Protected member functions

        //- Construct as copy (not implemented)
        basicPsiThermo(const basicPsiThermo&);


public:

    //- Runtime type information
    TypeName("basicPsiThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicPsiThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );

    // Constructors

        //- Construct from mesh
        basicPsiThermo(const fvMesh&);


    //- Selector
    static autoPtr<basicPsiThermo> New(const fvMesh&);

    //- Destructor
    virtual ~basicPsiThermo();


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3] - uses current value of pressure
            virtual tmp<volScalarField> rho() const
            {
                return p_*psi();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
