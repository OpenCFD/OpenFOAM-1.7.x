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
    Foam::duplicatePoints

Description
    Duplicate points

SourceFiles
    duplicatePoints.C

\*---------------------------------------------------------------------------*/

#ifndef duplicatePoints_H
#define duplicatePoints_H

#include "labelList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class edge;
class face;
class mapPolyMesh;
class localPointRegion;

/*---------------------------------------------------------------------------*\
                           Class duplicatePoints Declaration
\*---------------------------------------------------------------------------*/

class duplicatePoints
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Per regionSide.pointRegion() the points created.
        //  element 0 is the original point.
        //  Can have empty entries if a point only has one region.
        labelListList duplicates_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        duplicatePoints(const duplicatePoints&);

        //- Disallow default bitwise assignment
        void operator=(const duplicatePoints&);

public:

    //- Runtime type information
    ClassName("duplicatePoints");


    // Constructors

        //- Construct from mesh
        duplicatePoints(const polyMesh& mesh);


    // Member Functions

        // Access

            //- Per point in regionSide.pointRegions() the points created.
            const labelListList& duplicates() const
            {
                return duplicates_;
            }


        // Topology changes

            //- Play commands into polyTopoChange to duplicate points. Gets
            //  localPointRegion structure which is per non-manifold point
            //  the regions per point.
            void setRefinement
            (
                const localPointRegion& regionSide,
                polyTopoChange&
            );

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
