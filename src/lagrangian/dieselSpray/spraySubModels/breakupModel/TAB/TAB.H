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
    Foam::TAB

Description
    The TAB Method for Numerical Calculation of Spray Droplet Breakup.

    @verbatim
        O'Rourke, P.J. and Amsden, A.A.,
        "The TAB Method for Numerical Calculation of Spray Droplet Breakup,"
        1987 SAE International Fuels and Lubricants Meeting and Exposition,
        Toronto, Ontario, November 2-5, 1987,
        Los Alamos National Laboratory document LA-UR-87-2105;
        SAE Technical Paper Series, Paper 872089.
    @endverbatim

    This implementation follows the kiva version.

See Also
    The Enhanced %TAB model - ETAB

SourceFiles
    TABI.H
    TAB.C
    TABIO.C

\*---------------------------------------------------------------------------*/

#ifndef TAB_H
#define TAB_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TAB Declaration
\*---------------------------------------------------------------------------*/

class TAB
:
    public breakupModel
{
private:

    // Private data

        dictionary coeffsDict_;

    // inverse function approximation
    // of the Rossin-Rammler Distribution
    // used when calculating the droplet size after breakup
        scalar rrd_[100];

        // model constants

        scalar Cmu_;
        scalar Comega_;
        scalar WeCrit_;

public:


    //- Runtime type information
    TypeName("TAB");


    // Constructors

        //- Construct from components
        TAB
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~TAB();


    // Member Functions

        void breakupParcel
        (
            parcel& parcel,
            const scalar deltaT,
            const vector& Ug,
            const liquidMixture& fuels
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
