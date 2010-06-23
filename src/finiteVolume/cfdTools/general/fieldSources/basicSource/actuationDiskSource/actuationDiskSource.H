/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::actuationDiskSource

Description
    Actuation disk zone definition.
    Constant values for momentum source for actuation disk

    T = 2*rho*A*sqr(Uo)*a*(1-a)
    U1 = (1 -a)Uo
    where:
          A: disk area
          Uo: upstream velocity
          a: 1 - Cp/Ct
          U1: velocity at the disk

SourceFiles
    actuationDiskSource.C
    actuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef actuationDiskSource_H
#define actuationDiskSource_H

#include "IOdictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "DimensionedField.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "basicSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class actuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class actuationDiskSource
:
    public basicSource
{
    // Private data

        //- Cell zone ID
        label cellZoneID_;

        //- Sub dictionary with actuationDisk information
        const dictionary& dict_;

        //- Disk area normal
        vector diskDir_;

        //- Power coefficient
        scalar Cp_;

        //- Thrust coefficient
        scalar Ct_;

        //- Disk area
        scalar diskArea_;


    // Private Member Functions

        //- Check data
        void checkData();

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        actuationDiskSource(const actuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const actuationDiskSource&);


public:

    //- Runtime type information
    TypeName("actuationDiskSource");


    // Constructors

        //- Construct from components
        actuationDiskSource
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~actuationDiskSource()
    {}


    // Member Functions

        // Access

            //- cellZone number
            label zoneId() const
            {
                return cellZoneID_;
            }

            //- Return Cp
            scalar Cp() const
            {
                return Cp_;
            }

            //- Return Ct
            scalar Ct() const
            {
                return Ct_;
            }

            //- Normal disk direction
            const vector& diskDir() const
            {
                return diskDir_;
            }

            //- Disk area
            scalar diskArea() const
            {
                return diskArea_;
            }


        // Public Functions

            //-Source term to fvMatrix<vector>
            virtual void addSu(fvMatrix<vector>& UEqn);

            //-Source term to fvMatrix<scalar>
            virtual void addSu(fvMatrix<scalar>& UEqn){}

            //- Add all explicit source
            virtual void addExplicitSources(){}

            //- Add source to scalar field
            virtual void addSu(DimensionedField<scalar, volMesh>& field){}

            //- Add source to vector field
            virtual void addSu(DimensionedField<vector, volMesh>& field){}


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "actuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
