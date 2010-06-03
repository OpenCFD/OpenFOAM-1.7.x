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
    Foam::SHF

Description
    Secondary Breakup Model to take account of the different breakup regimes,
    bag, molutimode, shear....

    Accurate description in
    @verbatim
    R. Schmehl, G. Maier, S. Witting
    "CFD Analysis of Fuel Atomization, Secondary Droplet Breakup and Spray
    Dispersion in the Premix Duct of a LPP Combustor".
    Eight International Conference on Liquid Atomization and Spray Systems, 2000
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef SHF_H
#define SHF_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SHF Declaration
\*---------------------------------------------------------------------------*/

class SHF
:
    public breakupModel
{

private:

    // Private data

        dictionary coeffsDict_;

        // reference to gravity
        const vector& g_;

        Random& rndGen_;

        // model constants

        scalar weCorrCoeff_;

        scalar weBuCrit_;
        scalar weBuBag_;
        scalar weBuMM_;

        scalar ohnCoeffCrit_;
        scalar ohnCoeffBag_;
        scalar ohnCoeffMM_;

        scalar ohnExpCrit_;
        scalar ohnExpBag_;
        scalar ohnExpMM_;

        scalar cInit_;

        scalar c1_;
        scalar c2_;
        scalar c3_;

        scalar cExp1_;
        scalar cExp2_;
        scalar cExp3_;

        scalar weConst_;
        scalar weCrit1_;
        scalar weCrit2_;

        scalar coeffD_;
        scalar onExpD_;
        scalar weExpD_;

        scalar mu_;
        scalar sigma_;

        scalar d32Coeff_;
        scalar cDmaxBM_;
        scalar cDmaxS_;

        scalar corePerc_;

public:

    //- Runtime type information
    TypeName("SHF");


    // Constructors

        //- Construct from components
        SHF
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~SHF();


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
