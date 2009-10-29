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
    Foam::basicSensiblePsiThermo

Description
    Basic thermodynamic properties based on compressibility

SourceFiles
    basicSensiblePsiThermo.C
    newBasicSensiblePsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicSensiblePsiThermo_H
#define basicSensiblePsiThermo_H

#include "basicThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class basicSensiblePsiThermo Declaration
\*---------------------------------------------------------------------------*/

class basicSensiblePsiThermo
:
    public basicThermo
{

protected:

    // Protected member functions

        //- Construct as copy (not implemented)
        basicSensiblePsiThermo(const basicSensiblePsiThermo&);


public:

    //- Runtime type information
    TypeName("basicSensiblePsiThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicSensiblePsiThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );

    // Constructors

        //- Construct from mesh
        basicSensiblePsiThermo(const fvMesh&);


    //- Selector
    static autoPtr<basicSensiblePsiThermo> New(const fvMesh&);

    //- Destructor
    virtual ~basicSensiblePsiThermo();


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
