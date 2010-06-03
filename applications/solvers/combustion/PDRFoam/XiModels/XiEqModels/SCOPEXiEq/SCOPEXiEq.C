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

\*---------------------------------------------------------------------------*/

#include "SCOPEXiEq.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(SCOPEXiEq, 0);
    addToRunTimeSelectionTable(XiEqModel, SCOPEXiEq, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::SCOPEXiEq::SCOPEXiEq
(
    const dictionary& XiEqProperties,
    const hhuCombustionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, thermo, turbulence, Su),
    XiEqCoef(readScalar(XiEqModelCoeffs_.lookup("XiEqCoef"))),
    XiEqExp(readScalar(XiEqModelCoeffs_.lookup("XiEqExp"))),
    lCoef(readScalar(XiEqModelCoeffs_.lookup("lCoef"))),
    SuMin(0.01*Su.average()),
    MaModel
    (
        IOdictionary
        (
            IOobject
            (
                "combustionProperties",
                Su.mesh().time().constant(),
                Su.mesh(),
                IOobject::MUST_READ
            )
        ),
        thermo
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiEqModels::SCOPEXiEq::~SCOPEXiEq()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::SCOPEXiEq::XiEq() const
{
    const volScalarField& k = turbulence_.k();
    const volScalarField& epsilon = turbulence_.epsilon();

    volScalarField up = sqrt((2.0/3.0)*k);
    volScalarField l = (lCoef*sqrt(3.0/2.0))*up*k/epsilon;
    volScalarField Rl = up*l*thermo_.rhou()/thermo_.muu();

    volScalarField upBySu = up/(Su_ + SuMin);
    volScalarField K = 0.157*upBySu/sqrt(Rl);
    volScalarField Ma = MaModel.Ma();

    tmp<volScalarField> tXiEq
    (
        new volScalarField
        (
            IOobject
            (
                "XiEq",
                epsilon.time().timeName(),
                epsilon.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            epsilon.mesh(),
            dimensionedScalar("XiEq", dimless, 0.0)
        )
    );
    volScalarField& xieq = tXiEq();

    forAll(xieq, celli)
    {
        if (Ma[celli] > 0.01)
        {
            xieq[celli] =
                XiEqCoef*pow(K[celli]*Ma[celli], -XiEqExp)*upBySu[celli];
        }
    }

    forAll(xieq.boundaryField(), patchi)
    {
        scalarField& xieqp = xieq.boundaryField()[patchi];
        const scalarField& Kp = K.boundaryField()[patchi];
        const scalarField& Map = Ma.boundaryField()[patchi];
        const scalarField& upBySup = upBySu.boundaryField()[patchi];

        forAll(xieqp, facei)
        {
            if (Ma[facei] > 0.01)
            {
                xieqp[facei] =
                    XiEqCoef*pow(Kp[facei]*Map[facei], -XiEqExp)*upBySup[facei];
            }
        }
    }

    return tXiEq;
}


bool Foam::XiEqModels::SCOPEXiEq::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    XiEqModelCoeffs_.lookup("XiEqCoef") >> XiEqCoef;
    XiEqModelCoeffs_.lookup("XiEqExp") >> XiEqExp;
    XiEqModelCoeffs_.lookup("lCoef") >> lCoef;

    return true;
}


// ************************************************************************* //
