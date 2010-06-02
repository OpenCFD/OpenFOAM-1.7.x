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
    Foam::fv::cellMDLimitedGrad

Description
    cellMDLimitedGrad gradient scheme applied to a runTime selected base
    gradient scheme.

    The scalar limiter based on limiting the extrapolated face values
    between the maximum and minimum cell and cell neighbour values and is
    applied to the gradient in each face direction separately.

SourceFiles
    cellMDLimitedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef cellMDLimitedGrad_H
#define cellMDLimitedGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class cellMDLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellMDLimitedGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        tmp<fv::gradScheme<Type> > basicGradScheme_;

        //- Limiter coefficient
        const scalar k_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellMDLimitedGrad(const cellMDLimitedGrad&);

        //- Disallow default bitwise assignment
        void operator=(const cellMDLimitedGrad&);


public:

    //- RunTime type information
    TypeName("cellMDLimited");


    // Constructors

        //- Construct from mesh and schemeData
        cellMDLimitedGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            basicGradScheme_(fv::gradScheme<Type>::New(mesh, schemeData)),
            k_(readScalar(schemeData))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorIn
                (
                    "cellMDLimitedGrad(const fvMesh&, Istream& schemeData)",
                    schemeData
                )   << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        static inline void limitFace
        (
            typename outerProduct<vector, Type>::type& g,
            const Type& maxDelta,
            const Type& minDelta,
            const vector& dcf
        );

        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > grad
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
inline void cellMDLimitedGrad<scalar>::limitFace
(
    vector& g,
    const scalar& maxDelta,
    const scalar& minDelta,
    const vector& dcf
)
{
    scalar extrapolate = dcf & g;

    if (extrapolate > maxDelta)
    {
        g = g + dcf*(maxDelta - extrapolate)/magSqr(dcf);
    }
    else if (extrapolate < minDelta)
    {
        g = g + dcf*(minDelta - extrapolate)/magSqr(dcf);
    }
}


template<class Type>
inline void cellMDLimitedGrad<Type>::limitFace
(
    typename outerProduct<vector, Type>::type& g,
    const Type& maxDelta,
    const Type& minDelta,
    const vector& dcf
)
{
    for(direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        vector gi(g[cmpt], g[cmpt+3], g[cmpt+6]);
        cellMDLimitedGrad<scalar>::limitFace
        (
            gi,
            maxDelta.component(cmpt),
            minDelta.component(cmpt),
            dcf
        );
        g[cmpt] = gi.x();
        g[cmpt+3] = gi.y();
        g[cmpt+6] = gi.z();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
