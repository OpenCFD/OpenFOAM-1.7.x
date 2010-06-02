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
    Foam::XiEqModel

Description
    Base-class for all XiEq models used by the b-XiEq combustion model.
    The available models are :
                                \link basicXiSubXiEq.H \endlink
                                \link Gulder.H \endlink
                                \link instabilityXiEq.H \endlink
                                \link SCOPEBlendXiEq.H \endlink
                                \link SCOPEXiEq.H \endlink

SourceFiles
    XiEqModel.C

\*---------------------------------------------------------------------------*/

#ifndef XiEqModel_H
#define XiEqModel_H

#include "IOdictionary.H"
#include "hhuCombustionThermo.H"
#include "RASModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiEqModel Declaration
\*---------------------------------------------------------------------------*/

class XiEqModel
{

protected:

    // Protected data

        dictionary XiEqModelCoeffs_;

        const hhuCombustionThermo& thermo_;
        const compressible::RASModel& turbulence_;
        const volScalarField& Su_;


private:

    // Private Member Functions

        //- Disallow copy construct
        XiEqModel(const XiEqModel&);

        //- Disallow default bitwise assignment
        void operator=(const XiEqModel&);


public:

    //- Runtime type information
    TypeName("XiEqModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiEqModel,
            dictionary,
            (
                const dictionary& XiEqProperties,
                const hhuCombustionThermo& thermo,
                const compressible::RASModel& turbulence,
                const volScalarField& Su
            ),
            (
                XiEqProperties,
                thermo,
                turbulence,
                Su
            )
        );


    // Selectors

        //- Return a reference to the selected XiEq model
        static autoPtr<XiEqModel> New
        (
            const dictionary& XiEqProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Constructors

        //- Construct from components
        XiEqModel
        (
            const dictionary& XiEqProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Destructor

        virtual ~XiEqModel();


    // Member Functions

        //- Return the flame-wrinking XiEq
        virtual tmp<volScalarField> XiEq() const
        {
            return turbulence_.muEff();
        }

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
