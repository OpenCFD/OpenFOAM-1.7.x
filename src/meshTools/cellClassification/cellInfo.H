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
    Foam::cellInfo

Description
    Holds information regarding type of cell. Used in inside/outside 
    determination in cellClassification.

SourceFiles
    cellInfoI.H
    cellInfo.C

\*---------------------------------------------------------------------------*/

#ifndef cellInfo_H
#define cellInfo_H

#include "point.H"
#include "label.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class cellInfo Declaration
\*---------------------------------------------------------------------------*/

class cellInfo
{
    // Private data

        label type_;
        
    // Private Member Functions

        //- Update current cell/face type with neighbouring
        //  type. Return true if information needs to propagate,
        //  false otherwise.
        inline bool update
        (
            const cellInfo& w2,
            const label thisFaceI,
            const label thisCellI,
            const label neighbourFaceI,
            const label neighbourCellI
        );

public:

    // Constructors

        //- Construct null
        inline cellInfo();

        //- Construct from cType
        inline cellInfo(const label);

        //- Construct as copy
        inline cellInfo(const cellInfo&);


    // Member Functions

        // Access

            inline label type() const
            {
                return type_;
            }


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid() const;

            //- Check for identical geometrical data. Used for cyclics checking.
            inline bool sameGeometry
            (
                const polyMesh&,
                const cellInfo&,
                const scalar
            ) const;

            //- Convert any absolute coordinates into relative to (patch)face
            //  centre
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch& patch,
                const label patchFaceI,
                const point& faceCentre
            );

            //- Reverse of leaveDomain
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch& patch,
                const label patchFaceI,
                const point& faceCentre
            );

            //- Apply rotation matrix to any coordinates
            inline void transform
            (
                const polyMesh&,
                const tensor& rotTensor
            );

            //- Influence of neighbouring face.
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCellI,
                const label neighbourFaceI,
                const cellInfo& neighbourInfo,
                const scalar tol
            );

            //- Influence of neighbouring cell.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const label neighbourCellI,
                const cellInfo& neighbourInfo,
                const scalar tol
            );

            //- Influence of different value on same face.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const cellInfo& neighbourInfo,
                const scalar tol
            );

    // Member Operators

        //Note: Used to determine whether to call update.
        inline bool operator==(const cellInfo&) const;

        inline bool operator!=(const cellInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const cellInfo&);
        friend Istream& operator>>(Istream&, cellInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
