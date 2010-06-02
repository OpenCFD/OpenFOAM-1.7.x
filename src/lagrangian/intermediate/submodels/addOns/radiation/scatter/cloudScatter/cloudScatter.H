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
    Foam::radiation::cloudScatter

Description
    Cloud radiation scatter model

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef radiationCloudScatter_H
#define radiationCloudScatter_H

#include "scatterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                         Class cloudScatter Declaration
\*---------------------------------------------------------------------------*/

class cloudScatter
:
    public scatterModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Cloud names
        wordList cloudNames_;


public:

    //- Runtime type information
    TypeName("cloudScatter");


    // Constructors

        //- Construct from components
        cloudScatter
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Destructor
    virtual ~cloudScatter();


    // Member Operators

        // Access

            //- Return scatter coefficient
            tmp<volScalarField> sigmaEff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
