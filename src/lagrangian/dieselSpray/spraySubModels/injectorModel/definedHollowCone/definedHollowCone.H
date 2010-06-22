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
    Foam::definedHollowConeInjector

Description
    An injector model for hollow/solid cone injection

    - user defined:
      - inner/outer cone angle vs. time
      - droplet PDF

\*---------------------------------------------------------------------------*/

#ifndef definedHollowConeInjector_H
#define definedHollowConeInjector_H

#include "injectorModel.H"
#include "vector.H"
#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class definedHollowConeInjector Declaration
\*---------------------------------------------------------------------------*/

class definedHollowConeInjector
:
    public injectorModel
{

private:

    // Private data

        typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

        dictionary definedHollowConeDict_;
        autoPtr<pdfs::pdf> dropletPDF_;

        // inner and outer cone angle time histories
        // 2 column vectors = (time, coneAngle)
        List<pair> innerConeAngle_;
        List<pair> outerConeAngle_;

public:

    //- Runtime type information
    TypeName("definedHollowConeInjector");


    // Constructors

        //- Construct from components
        definedHollowConeInjector
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~definedHollowConeInjector();


    // Member Functions

        //- Return the injected droplet diameter
        scalar d0(const label injector, const scalar time) const;

        //- Return the spray angle of the injector
        vector direction
        (
            const label injector,
            const label hole,
            const scalar time,
            const scalar d
        ) const;

        scalar velocity
        (
            const label i,
            const scalar time
        ) const;

        scalar averageVelocity
        (
            const label i
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
