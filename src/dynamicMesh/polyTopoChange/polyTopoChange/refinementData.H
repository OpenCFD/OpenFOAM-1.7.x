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
    Foam::refinementData

Description
    Transfers refinement levels such that slow transition between levels is
    maintained. Used in meshWave.

SourceFiles
    refinementDataI.H
    refinementData.C

\*---------------------------------------------------------------------------*/

#ifndef refinementData_H
#define refinementData_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class refinementData Declaration
\*---------------------------------------------------------------------------*/

class refinementData
{
    // Private data

        //- Count which triggers refinement
        label refinementCount_;

        //- Refinement level
        label count_;

public:

    // Constructors

        //- Construct null
        inline refinementData();

        //- Construct from count
        inline refinementData(const label refinementCount, const label count);


    // Member Functions

        // Access

            inline label refinementCount() const
            {
                return refinementCount_;
            }

            inline label& refinementCount()
            {
                return refinementCount_;
            }

            inline label count() const
            {
                return count_;
            }

            inline label& count()
            {
                return count_;
            }

            inline bool isRefined() const
            {
                return count_ >= refinementCount_;
            }



        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid() const;

            //- Check for identical geometrical data. Used for cyclics checking.
            inline bool sameGeometry
            (
                const polyMesh&,
                const refinementData&,
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
                const refinementData& neighbourInfo,
                const scalar tol
            );

            //- Influence of neighbouring cell.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const label neighbourCellI,
                const refinementData& neighbourInfo,
                const scalar tol
            );

            //- Influence of different value on same face.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const refinementData& neighbourInfo,
                const scalar tol
            );

    // Member Operators

        // Needed for List IO
        inline bool operator==(const refinementData&) const;

        inline bool operator!=(const refinementData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const refinementData&);
        friend Istream& operator>>(Istream&, refinementData&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "refinementDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
