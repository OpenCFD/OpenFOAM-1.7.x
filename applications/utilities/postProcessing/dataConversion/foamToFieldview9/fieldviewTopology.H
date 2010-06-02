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
    Foam::fieldviewTopology

Description

SourceFiles
    fieldviewTopology.C

\*---------------------------------------------------------------------------*/

#ifndef fieldviewTopology_H
#define fieldviewTopology_H

#include "labelList.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class cellShape;

/*---------------------------------------------------------------------------*\
                           Class fieldviewTopology Declaration
\*---------------------------------------------------------------------------*/

class fieldviewTopology
{
    // Private data

        //- Hexes in fieldview format
        labelList hexLabels_;

        labelList prismLabels_;

        labelList pyrLabels_;

        labelList tetLabels_;

        //- Number of non-hex/prism/pyr/tet labels
        label nPoly_;


        //
        // Patches
        //

        //- Quad and tri patch faces in fv format
        labelListList quadFaceLabels_;

        //- Number of polyhedral faces per patch
        labelList nPolyFaces_;


    // Private Member Functions

        static labelList calcFaceAddressing
        (
            const faceList& allFaces,   // faces given faceLabels
            const cellShape& shape,
            const labelList& faces,     // faceLabels for given cell
            const label cellI
        );


        //- Disallow default bitwise copy construct
        fieldviewTopology(const fieldviewTopology&);

        //- Disallow default bitwise assignment
        void operator=(const fieldviewTopology&);


public:

    // Constructors

        //- Construct from components
        fieldviewTopology(const polyMesh& mesh, const bool setWallInfo);


    // Member Functions

        // Access

            const labelList& hexLabels() const
            {
                return hexLabels_;
            }

            const labelList& prismLabels() const
            {
                return prismLabels_;
            }

            const labelList& pyrLabels() const
            {
                return pyrLabels_;
            }

            const labelList& tetLabels() const
            {
                return tetLabels_;
            }

            label nHex() const
            {
                return hexLabels().size()/9;
            }

            label nPrism() const
            {
                return prismLabels().size()/7;
            }

            label nPyr() const
            {
                return pyrLabels().size()/6;
            }

            label nTet() const
            {
                return tetLabels().size()/5;
            }

            label nPoly() const
            {
                return nPoly_;
            }

            const labelListList& quadFaceLabels() const
            {
                return quadFaceLabels_;
            }

            const labelList& nPolyFaces() const
            {
                return nPolyFaces_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
