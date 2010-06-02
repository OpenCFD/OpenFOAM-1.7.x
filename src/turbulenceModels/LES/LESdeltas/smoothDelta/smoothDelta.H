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
    Foam::smoothDelta

Description
    Smoothed delta which takes a given simple geometric delta and applies
    smoothing to it such that the ratio of deltas between two cells is no
    larger than a specified amount, typically 1.15.

SourceFiles
    smoothDelta.C

\*---------------------------------------------------------------------------*/

#ifndef smoothDelta_H
#define smoothDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class smoothDelta Declaration
\*---------------------------------------------------------------------------*/

class smoothDelta
:
    public LESdelta
{
    // Private data

        autoPtr<LESdelta> geometricDelta_;
        scalar maxDeltaRatio_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        smoothDelta(const smoothDelta&);
        void operator=(const smoothDelta&);

        // Calculate the delta values
        void calcDelta();

        //- Private member class used by mesh-wave to propagate the delta-ratio
        class deltaData
        {
            scalar delta_;

            // Private Member Functions

                //- Update. Gets information from neighbouring face/cell and
                //  uses this to update itself (if nessecary) and return true.
                inline bool update
                (
                    const deltaData& w2,
                    const scalar scale,
                    const scalar tol
                );


        public:

            // Static data members

                //- delta fraction
                static scalar maxDeltaRatio;


                // Constructors

                //- Construct null
                inline deltaData();

                //- Construct from origin, yStar, distance
                inline deltaData(const scalar delta);


            // Member Functions

                // Access

                scalar delta() const
                {
                    return delta_;
                }


                // Needed by FaceCellWave

                    //- Check whether origin has been changed at all or
                    //  still contains original (invalid) value.
                    inline bool valid() const;

                    //- Check for identical geometrical data.
                    //  Used for cyclics checking.
                    inline bool sameGeometry
                    (
                        const polyMesh&,
                        const deltaData&,
                        const scalar
                    ) const;

                    //- Convert any absolute coordinates into relative to
                    //  (patch)face centre
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
                        const deltaData& neighbourInfo,
                        const scalar tol
                    );

                    //- Influence of neighbouring cell.
                    inline bool updateFace
                    (
                        const polyMesh&,
                        const label thisFaceI,
                        const label neighbourCellI,
                        const deltaData& neighbourInfo,
                        const scalar tol
                    );

                    //- Influence of different value on same face.
                    inline bool updateFace
                    (
                        const polyMesh&,
                        const label thisFaceI,
                        const deltaData& neighbourInfo,
                        const scalar tol
                    );

                // Member Operators

                    // Needed for List IO
                    inline bool operator==(const deltaData&) const;

                    inline bool operator!=(const deltaData&) const;

                // IOstream Operators

                    friend Ostream& operator<<
                    (
                        Ostream& os,
                        const deltaData& wDist
                    )
                    {
                        return os << wDist.delta_;
                    }

                    friend Istream& operator>>(Istream& is, deltaData& wDist)
                    {
                        return is >> wDist.delta_;
                    }
        };


        void setChangedFaces
        (
            const polyMesh& mesh,
            const volScalarField& delta,
            DynamicList<label>& changedFaces,
            DynamicList<deltaData>& changedFacesInfo
        );


public:

    //- Runtime type information
    TypeName("smooth");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        smoothDelta
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        );


    //- Destructor
    virtual ~smoothDelta()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "smoothDeltaDeltaDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
