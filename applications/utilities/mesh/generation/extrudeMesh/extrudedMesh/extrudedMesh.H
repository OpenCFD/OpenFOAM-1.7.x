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
    Foam::extrudedMesh

Description

SourceFiles
    extrudedMesh.C

\*---------------------------------------------------------------------------*/

#ifndef extrudedMesh_H
#define extrudedMesh_H

#include "polyMesh.H"
#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class extrudedMesh Declaration
\*---------------------------------------------------------------------------*/

class extrudedMesh
:
    public polyMesh
{
    // Private data

        const extrudeModel& model_;


    // Private Member Functions

        //- Do edge and face use points in same order?
        static bool sameOrder(const face&, const edge&);

        //- Construct and return the extruded mesh points
        template<class Face, template<class> class FaceList, class PointField>
        Xfer<pointField> extrudedPoints
        (
            const PrimitivePatch<Face, FaceList, PointField>& extrudePatch,
            const extrudeModel&
        );

        //- Construct and return the extruded mesh faces
        template<class Face, template<class> class FaceList, class PointField>
        Xfer<faceList> extrudedFaces
        (
            const PrimitivePatch<Face, FaceList, PointField>& extrudePatch,
            const extrudeModel&
        );

        //- Construct and return the extruded mesh cells
        template<class Face, template<class> class FaceList, class PointField>
        Xfer<cellList> extrudedCells
        (
            const PrimitivePatch<Face, FaceList, PointField>& extrudePatch,
            const extrudeModel&
        );


        //- Disallow default bitwise copy construct
        extrudedMesh(const extrudedMesh&);

        //- Disallow default bitwise assignment
        void operator=(const extrudedMesh&);


public:

    // Constructors

        //- Construct from the primitivePatch to extrude
        template<class Face, template<class> class FaceList, class PointField>
        extrudedMesh
        (
            const IOobject&,
            const PrimitivePatch<Face, FaceList, PointField>& extrudePatch,
            const extrudeModel&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extrudedMesh.C"
#else
#   ifdef xlC
#   pragma implementation("extrudedMesh.C")
#   endif
#endif

#endif

// ************************************************************************* //

