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
    Foam::refinementDistanceData

Description
    Transfers refinement levels such that slow transition between levels is
    maintained. Used in FaceCellWave.

SourceFiles
    refinementDistanceDataI.H
    refinementDistanceData.C

\*---------------------------------------------------------------------------*/

#ifndef refinementDistanceData_H
#define refinementDistanceData_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class refinementDistanceData Declaration
\*---------------------------------------------------------------------------*/

class refinementDistanceData
{

    // Private data

        //- unrefined (level0) buffer size (nBufferLayers*level0Size)
        scalar level0Size_;

        //- nearest point with highest level
        point origin_;
        label originLevel_;


    // Private Member Functions

        //- Updates with neighbouring data. Returns true if something changed.
        inline bool update
        (
            const point&,
            const refinementDistanceData& neighbourInfo,
            const scalar tol
        );

public:

    // Constructors

        //- Construct null
        inline refinementDistanceData();

        //- Construct from count
        inline refinementDistanceData
        (
            const scalar level0Size,
            const point& origin,
            const label level
        );


    // Member Functions

        // Access

            inline scalar level0Size() const
            {
                return level0Size_;
            }

            inline scalar& level0Size()
            {
                return level0Size_;
            }

            inline const point& origin() const
            {
                return origin_;
            }

            inline point& origin()
            {
                return origin_;
            }

            inline label originLevel() const
            {
                return originLevel_;
            }

            inline label& originLevel()
            {
                return originLevel_;
            }


        // Other

            //- Calculates the wanted level at a given point. Walks out from
            //  the origin.
            inline label wantedLevel(const point& pt) const;


        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid() const;

            //- Check for identical geometrical data. Used for cyclics checking.
            inline bool sameGeometry
            (
                const polyMesh&,
                const refinementDistanceData&,
                const scalar
            ) const;

            //- Convert any absolute coordinates into relative to (patch)face
            //  centre
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFaceI,
                const point& faceCentre
            );

            //- Reverse of leaveDomain
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFaceI,
                const point& faceCentre
            );

            //- Apply rotation matrix to any coordinates
            inline void transform
            (
                const polyMesh&,
                const tensor&
            );

            //- Influence of neighbouring face.
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCellI,
                const label neighbourFaceI,
                const refinementDistanceData& neighbourInfo,
                const scalar tol
            );

            //- Influence of neighbouring cell.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const label neighbourCellI,
                const refinementDistanceData& neighbourInfo,
                const scalar tol
            );

            //- Influence of different value on same face.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const refinementDistanceData& neighbourInfo,
                const scalar tol
            );

    // Member Operators

        // Needed for List IO
        inline bool operator==(const refinementDistanceData&) const;

        inline bool operator!=(const refinementDistanceData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const refinementDistanceData&);
        friend Istream& operator>>(Istream&, refinementDistanceData&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "refinementDistanceDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
