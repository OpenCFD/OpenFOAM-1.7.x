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
    Foam::timeActivatedExplicitMulticomponentPointSourceNew

Description
    Provides a mechanism to introduce point sources to a set of carrier fields.
    Carrier fields are supplied on consruction, and interrogated to provide the
    field indices of the sources.

    Source values are assumed to be in <quantity>/s, and divided through by the
    cell volumes before being returned, e.g. for a kg/s source

    Properties are described in a <name>Properties dictionary

    active          true;       // are sources active (true/false)

    pointSources
    (
        source1                 // source name
        {
            timeStart           0.0;
            duration            1.0;
            location            (0 0 0);
            fieldData
            (
                (H2O            0.1)               // kg/s
                (O2             0.05)              // kg/s
            );
        }
        source2                 // source name
        {
            timeStart           0.5;
            duration            2.0;
            location            (1 1 1);
            fieldData
            (
                (NO             0.1)               // kg/s
                (CO2            0.05)              // kg/s
                (H2             0.001)             // kg/s
            );
        }
    );


SourceFiles
    timeActivatedExplicitMulticomponentPointSourceNew.C

\*---------------------------------------------------------------------------*/

#ifndef timeActivatedExplicitMulticomponentPointSource_H
#define timeActivatedExplicitMulticomponentPointSource_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "pointSourceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class timeActivatedExplicitMulitcomponentPointSource Declaration
\*---------------------------------------------------------------------------*/

class timeActivatedExplicitMulticomponentPointSource
:
    public IOdictionary
{
protected:

    // Protected data

        //- Name of the source
        word name_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to time database
        const Time& runTime_;

        //- Source dimensions
        const dimensionSet& dimensions_;

        //- Reference to the multicomponent carrier fields
        const PtrList<volScalarField>& carrierFields_;

        //- Active flag
        Switch active_;

        //- List of point source properties
        List<pointSourceProperties> pointSources_;

        //- List of cell owners for point source locations
        List<label> cellOwners_;

        //- List of field ids for each source
        List<labelList> fieldIds_;


    // Protected Member Functions

        //- Return the id of field given its name
        label carrierFieldId(const word& fieldName);

        //- Update the addressing between source and carrier fields
        void updateAddressing();

        //- Disallow default bitwise copy construct
        timeActivatedExplicitMulticomponentPointSource
        (
            const timeActivatedExplicitMulticomponentPointSource&
        );

        //- Disallow default bitwise assignment
        void operator=(const timeActivatedExplicitMulticomponentPointSource&);


public:

    // Constructors

        //- Construct from components
        timeActivatedExplicitMulticomponentPointSource
        (
            const word&,
            const fvMesh&,
            const PtrList<volScalarField>&,
            const dimensionSet&
        );


    // Member Functions

        // Access

            //- Return a tmp field of the source for field fieldI
            virtual tmp<DimensionedField<scalar, volMesh> > Su
            (
                const label fieldI
            );

            //- Return a tmp field of the total source
            virtual tmp<DimensionedField<scalar, volMesh> > Su();


        //- Read properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
