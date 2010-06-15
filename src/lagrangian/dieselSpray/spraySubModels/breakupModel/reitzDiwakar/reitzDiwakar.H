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
    Foam::reitzDiwakar

Description
    secondary breakup model for the spray

\*---------------------------------------------------------------------------*/

#ifndef reitzDiwakar_H
#define reitzDiwakar_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reitzDiwakar Declaration
\*---------------------------------------------------------------------------*/

class reitzDiwakar
:
    public breakupModel
{
private:

    // Private data

        dictionary coeffsDict_;
        scalar Cbag_;
        scalar Cb_;
        scalar Cstrip_;
        scalar Cs_;

public:

    //- Runtime type information
    TypeName("ReitzDiwakar");


    // Constructors

        //- Construct from components
        reitzDiwakar
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~reitzDiwakar();


    // Member Operators

        void breakupParcel
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
