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
    Foam::porousZones

Description
    A centralized porousZone collection.

    Container class for a set of porousZones with the porousZone member
    functions implemented to loop over the functions for each porousZone.

    The input file @c constant/porousZone is implemented as
    IOPtrList\<porousZone\> and contains the following type of data:

    @verbatim
    1
    (
    cat1
    {
        coordinateSystem    system_10;
        porosity    0.781;
        Darcy
        {
            d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
            f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
        }
    }
    )
    @endverbatim

SourceFiles
    porousZones.C

\*---------------------------------------------------------------------------*/

#ifndef porousZones_H
#define porousZones_H

#include "porousZone.H"
#include "IOPtrList.H"

#include "volFieldsFwd.H"
#include "fvMatrix.H"
#include "geometricOneField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class porousZones Declaration
\*---------------------------------------------------------------------------*/

class porousZones
:
    public IOPtrList<porousZone>
{
    // Private data

        //- Reference to the finite volume mesh this zone is part of
        const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        porousZones(const porousZones&);

        //- Disallow default bitwise assignment
        void operator=(const porousZones&);


        //- modify time derivative elements
        template<class Type>
        void modifyDdt(fvMatrix<Type>&) const;

public:

    // Constructors

        //- Construct from fvMesh
        //  with automatically constructed coordinate systems list
        porousZones(const fvMesh&);


    // Member Functions

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type> > ddt
        (
            GeometricField<Type, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type> > ddt
        (
            const geometricOneField&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type> > ddt
        (
            const dimensionedScalar&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type> > ddt
        (
            const volScalarField&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        void addResistance(fvVectorMatrix& UEqn) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;

        //- read modified data
        virtual bool readData(Istream&);

        //- write data
        bool writeData(Ostream&, bool subDict = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "porousZonesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
