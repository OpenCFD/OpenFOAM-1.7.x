/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "SpalartAllmarasIDDES.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmarasIDDES, 0);
addToRunTimeSelectionTable(LESModel, SpalartAllmarasIDDES, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmarasIDDES::alpha() const
{
    return max
    (
        0.25 - y_/dimensionedScalar("hMax", dimLength, max(cmptMax(delta()))),
        scalar(-5)
    );
}


tmp<volScalarField> SpalartAllmarasIDDES::ft
(
    const volScalarField& S
) const
{
    return tanh(pow3(sqr(ct_)*rd(nuSgs_, S)));
}


tmp<volScalarField> SpalartAllmarasIDDES::fl
(
    const volScalarField& S
) const
{
    return tanh(pow(sqr(cl_)*rd(nu(), S), 10));
}


tmp<volScalarField> SpalartAllmarasIDDES::rd
(
    const volScalarField& visc,
    const volScalarField& S
) const
{
    return min
    (
        visc
       /(
           max
           (
               S,
               dimensionedScalar("SMALL", S.dimensions(), SMALL)
           )*sqr(kappa_*y_)
         + dimensionedScalar
           (
               "ROOTVSMALL",
               dimensionSet(0, 2 , -1, 0, 0),
               ROOTVSMALL
           )
       ),
       scalar(10)
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmarasIDDES::fd(const volScalarField& S) const
{
    return 1 - tanh(pow3(8*rd(nuEff(), S)));
}


tmp<volScalarField> SpalartAllmarasIDDES::dTilda(const volScalarField& S) const
{
    volScalarField alpha = this->alpha();
    volScalarField expTerm = exp(sqr(alpha));

    volScalarField fHill =
        2*(pos(alpha)*pow(expTerm, -11.09) + neg(alpha)*pow(expTerm, -9.0));

    volScalarField fStep = min(2*pow(expTerm, -9.0), scalar(1));
    volScalarField fHyb = max(1 - fd(S), fStep);
    volScalarField fAmp = 1 - max(ft(S), fl(S));
    volScalarField fRestore = max(fHill - 1, scalar(0))*fAmp;

    // IGNORING ft2 terms
    volScalarField Psi = sqrt
    (
        min
        (
            scalar(100),
            (1 - Cb1_/(Cw1_*sqr(kappa_)*fwStar_)*fv2())/max(SMALL, fv1())
        )
    );

    return max
    (
        dimensionedScalar("SMALL", dimLength, SMALL),
        fHyb*(1 + fRestore*Psi)*y_
      + (1 - fHyb)*CDES_*Psi*delta()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmarasIDDES::SpalartAllmarasIDDES
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    SpalartAllmaras(U, phi, transport, typeName),

    fwStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fwStar",
            coeffDict_,
            0.424
        )
    ),
    cl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cl",
            coeffDict_,
            3.55
        )
    ),
    ct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ct",
            coeffDict_,
            1.63
        )
    )

{}


bool SpalartAllmarasIDDES::read()
{
    if (SpalartAllmaras::read())
    {
        fwStar_.readIfPresent(coeffDict());
        cl_.readIfPresent(coeffDict());
        ct_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
