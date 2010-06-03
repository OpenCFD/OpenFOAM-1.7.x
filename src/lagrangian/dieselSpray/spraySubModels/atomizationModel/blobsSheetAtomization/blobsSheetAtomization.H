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
    Foam::blobsSheetAtomization

Description
    Primary Breakup Model for pressure swirl atomizers.

    Accurate description in
    @verbatim
    Z. Han, S. Parrish, P.V. Farrell, R.D. Reitz
    "Modeling Atomization Processes Of Pressure Swirl Hollow-Cone Fuel Sprays"
    Atomization and Sprays, vol. 7, pp. 663-684, 1997

    and

    L. Allocca, G. Bella, A. De Vita, L. Di Angelo
    "Experimental Validation of a GDI Spray Model"
    SAE Technical Paper Series, 2002-01-1137
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef blobsSheetAtomization_H
#define blobsSheetAtomization_H

#include "atomizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blobsSheetAtomization Declaration
\*---------------------------------------------------------------------------*/

class blobsSheetAtomization
:
    public atomizationModel
{
private:

    // Private data

        dictionary coeffsDict_;
        scalar B_;
        scalar angle_;
        Random& rndGen_;

public:

    //- Runtime type information
    TypeName("blobsSheetAtomization");


    // Constructors

        //- Construct from components
        blobsSheetAtomization
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~blobsSheetAtomization();


    // Member Operators

        void atomizeParcel
        (
            parcel& parcel,
            const scalar deltaT,
            const vector& vel,
            const liquidMixture& fuels
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
