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
    Foam::wallModel

Description
    A wall interaction model

SourceFiles
    wallModelI.H
    wallModel.C
    wallModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef wallModel_H
#define wallModel_H

#include "spray.H"
#include "IOdictionary.H"
#include "fvCFD.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wallModel Declaration
\*---------------------------------------------------------------------------*/

class wallModel
{

protected:

    // Protected data

        const dictionary& dict_;
        const fvMesh& mesh_;
        spray& spray_;


public:

    //- Runtime type information
        TypeName("wallModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            wallModel,
            dictionary,
            (
                const dictionary& dict,
                const volVectorField& U,
                spray& sm
            ),
            (dict, U, sm)
        );


    // Constructors

        //- Construct from components
        wallModel
        (
            const dictionary& dict,
            const volVectorField& U,
            spray& sm
        );


    // Destructor

        virtual ~wallModel();


    // Selector

        static autoPtr<wallModel> New
        (
            const dictionary& dict,
            const volVectorField& U,
            spray& sm
        );


    // Member Functions

        //- Return true if parcel is to be kept, false if it is to be removed
        virtual bool wallTreatment
        (
            parcel& parcel,
            const label facei
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
