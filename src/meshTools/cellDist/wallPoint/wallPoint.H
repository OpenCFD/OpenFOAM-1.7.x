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
    Foam::wallPoint

Description
    Holds information regarding nearest wall point. Used in wall distance
    calculation.

SourceFiles
    wallPointI.H
    wallPoint.C

\*---------------------------------------------------------------------------*/

#ifndef wallPoint_H
#define wallPoint_H

#include "point.H"
#include "label.H"
#include "scalar.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class wallPoint;

// Forward declaration of friend functions and operators
Ostream& operator<<(Ostream&, const wallPoint&);
Istream& operator>>(Istream&, wallPoint&);


/*---------------------------------------------------------------------------*\
                           Class wallPoint Declaration
\*---------------------------------------------------------------------------*/

class wallPoint
{
    // Private data

        //- position of nearest wall center
        point origin_;

        //- normal distance (squared) from cellcenter to origin
        scalar distSqr_;

    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        inline bool update
        (
            const point&,
            const wallPoint& w2,
            const scalar tol
        );


public:

    // Static data members

        //- initial point far away.
        static point greatPoint;


    // Constructors

        //- Construct null
        inline wallPoint();

        //- Construct from origin, distance
        inline wallPoint
        (
            const point& origin,
            const scalar distSqr
        );

        //- Construct as copy
        inline wallPoint
        (
            const wallPoint&
        );


    // Member Functions

        // Access

            inline const point& origin() const;

            inline point& origin();

            inline scalar distSqr() const;

            inline scalar& distSqr();


        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid() const;

            //- Check for identical geometrical data. Used for cyclics checking.
            inline bool sameGeometry
            (
                const polyMesh&,
                const wallPoint&,
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
                const wallPoint& neighbourInfo,
                const scalar tol
            );

            //- Influence of neighbouring cell.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const label neighbourCellI,
                const wallPoint& neighbourInfo,
                const scalar tol
            );

            //- Influence of different value on same face.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const wallPoint& neighbourInfo,
                const scalar tol
            );


    // Member Operators

        // Needed for List IO
        inline bool operator==(const wallPoint&) const;
        inline bool operator!=(const wallPoint&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const wallPoint&);
        friend Istream& operator>>(Istream&, wallPoint&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
