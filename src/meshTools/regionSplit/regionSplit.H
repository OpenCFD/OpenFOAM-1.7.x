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
    Foam::regionSplit

Description
    This class separates the mesh into distinct unconnected regions,
    each of which is then given a label.

SourceFiles
    regionSplit.C

\*---------------------------------------------------------------------------*/

#ifndef regionSplit_H
#define regionSplit_H

#include "polyMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regionSplit Declaration
\*---------------------------------------------------------------------------*/

class regionSplit
:
    public labelList
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Number of regions
        label nRegions_;


    // Private Member Functions

        //- Transfer faceRegion data from one face to the other (or vice versa)
        void transferCoupledFaceRegion
        (
            const label faceI,
            const label otherFaceI,

            labelList& faceRegion,
            DynamicList<label>& newChangedFaces
        ) const;

        //- Given a seed cell label, fill cellRegion/faceRegion with markValue
        //  for contiguous region around it
        void fillSeedMask
        (
            const List<labelPair>& explicitConnections,
            labelList& cellRegion,
            labelList& faceRegion,
            const label seedCellID,
            const label markValue
        ) const;

        //- Calculate region split. Return number of regions.
        label calcRegionSplit
        (
            const boolList& blockedFace,
            const List<labelPair>& explicitConnections,
            labelList& cellRegion
        ) const;

public:

    //- Runtime type information
    ClassName("regionSplit");

    // Constructors

        //- Construct from mesh
        regionSplit(const polyMesh&);

        //- Construct from mesh and whether face is blocked
        //  NOTE: blockedFace has to be consistent across coupled faces!
        regionSplit(const polyMesh&, const boolList& blockedFace);

        //- Construct from mesh and whether face is blocked. Additional explicit
        //  connections between normal boundary faces.
        //  NOTE: blockedFace has to be consistent across coupled faces!
        regionSplit
        (
            const polyMesh&,
            const boolList& blockedFace,
            const List<labelPair>&
        );

    // Member Functions

        //- Return number of regions
        label nRegions() const
        {
            return nRegions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
