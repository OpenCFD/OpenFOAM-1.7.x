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

Namespace
    Foam::meshWriters

Description
    A namespace for holding various types of mesh writers.


Class
    Foam::meshWriter

Description
    write OpenFOAM meshes and/or results to another CFD format
    - currently just STAR-CD

@par Files

    "constant/boundaryRegion" is an IOMap<dictionary> that contains
    the boundary type and names. eg,
    @verbatim
        (
            0
            {
                BoundaryType    wall;
                Label           Default_Boundary_Region;
            }

            1
            {
                BoundaryType    inlet;
                Label           inlet_1;
            }

            ...

            4
            {
                BoundaryType    pressure;
                Label           outlet;
            }
        )
    @endverbatim


SourceFiles
    meshWriterI.H
    meshWriter.C
    meshWriterIO.C

\*---------------------------------------------------------------------------*/

#ifndef meshWriter_H
#define meshWriter_H

#include "polyMesh.H"
#include "boundaryRegion.H"
#include "cellTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class meshWriter Declaration
\*---------------------------------------------------------------------------*/

class meshWriter
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshWriter(const meshWriter&);

        //- Disallow default bitwise assignment
        void operator=(const meshWriter&);


protected:

    // Protected data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Scaling factor for points (eg, [m] -> [mm])
        scalar scaleFactor_;

        //- Write bnd file
        bool writeBoundary_;

        //- boundaryRegion persistent data saved as a dictionary
        boundaryRegion boundaryRegion_;

        //- cellTable persistent data saved as a dictionary
        cellTable cellTable_;

        //- cellTable IDs for each cell
        labelList cellTableId_;

        //- Pointers to cell shape models
        static const cellModel* unknownModel;
        static const cellModel* tetModel;
        static const cellModel* pyrModel;
        static const cellModel* prismModel;
        static const cellModel* hexModel;


public:

    // Static data members

        static string defaultMeshName;
        static string defaultSurfaceName;


    // Constructors

        //- Ccreate a writer obejct
        meshWriter
        (
            const polyMesh&,
            const scalar scaleFactor = 1.0
        );


    //- Destructor
    virtual ~meshWriter();


    // Member Functions

        // Edit

            //- Set points scaling
            void scaleFactor(const scalar scaling)
            {
                scaleFactor_ = scaling;
            }

            //- Suppress writing bnd file
            void noBoundary()
            {
                writeBoundary_ = false;
            }


        // Write

            //- Write volume mesh
            //  subclass must to supply this method
            virtual bool write
            (
                const fileName& timeName = fileName::null
            ) const = 0;

            //- Write surface mesh with optional triangulation
            //  subclass could supply this information
            virtual bool writeSurface
            (
                const fileName& timeName = fileName::null,
                const bool& triangulate = false
            ) const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
