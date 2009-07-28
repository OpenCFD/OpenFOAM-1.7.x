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
    Foam::foamChemistryReader

Description
    Chemistry reader for OpenFOAM format

SourceFiles
    foamChemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef foamChemistryReader_H
#define foamChemistryReader_H

#include "chemistryReader.H"
#include "fileName.H"
#include "typeInfo.H"
#include "HashPtrTable.H"
#include "SLPtrList.H"
#include "labelList.H"
#include "speciesTable.H"
#include "atomicWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class foamChemistry Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class foamChemistryReader
:
    public chemistryReader<ThermoType>
{
        //- Table of the thermodynamic data given in the foamChemistry file
        HashPtrTable<ThermoType> speciesThermo_;

        //- Table of species
        speciesTable speciesTable_;

        //- List of the reactions
        SLPtrList<Reaction<ThermoType> > reactions_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        foamChemistryReader(const foamChemistryReader&);

        //- Disallow default bitwise assignment
        void operator=(const foamChemistryReader&);


public:

    //- Runtime type information
    TypeName("foamChemistryReader");


    // Constructors

        //- Construct from foamChemistry and thermodynamics file names
        foamChemistryReader
        (
            const fileName& reactionsFileName,
            const fileName& thermoFileName
        );

        //- Construct by getting the foamChemistry and thermodynamics file names
        //  from dictionary
        foamChemistryReader(const dictionary& thermoDict);


    //- Destructor
    virtual ~foamChemistryReader()
    {}


    // Member functions

        //- Table of species
        const speciesTable& species() const
        {
            return speciesTable_;
        }

        //- Table of the thermodynamic data given in the foamChemistry file
        const HashPtrTable<ThermoType>& speciesThermo() const
        {
            return speciesThermo_;
        }

        //- List of the reactions
        const SLPtrList<Reaction<ThermoType> >& reactions() const
        {
            return reactions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "foamChemistryReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
