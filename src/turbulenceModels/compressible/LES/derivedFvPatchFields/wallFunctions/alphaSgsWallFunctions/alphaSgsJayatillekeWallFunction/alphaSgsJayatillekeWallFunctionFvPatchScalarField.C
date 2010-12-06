/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "alphaSgsJayatillekeWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphaSgsJayatillekeWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphaSgsJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphaSgsJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void alphaSgsJayatillekeWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "alphaSgsJayatillekeWallFunctionFvPatchScalarField::checkType()"
        )
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


scalar alphaSgsJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar alphaSgsJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11.0;

    for (int i=0; i<maxIters_; i++)
    {
        scalar f = ypt - (log(E_*ypt)/kappa_ + P)/Prat;
        scalar df = 1.0 - 1.0/(ypt*kappa_*Prat);
        scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance_)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
     }

    return ypt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85),
    kappa_(0.41),
    E_(9.8),
    hsName_("hs")
{
    checkType();
}


alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const alphaSgsJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    hsName_(ptf.hsName_)
{}


alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    hsName_(dict.lookupOrDefault<word>("hs", "hs"))
{
    checkType();
}


alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const alphaSgsJayatillekeWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_),
    hsName_(awfpsf.hsName_)
{
    checkType();
}


alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const alphaSgsJayatillekeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_),
    hsName_(awfpsf.hsName_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphaSgsJayatillekeWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const LESModel& lesModel = db().lookupObject<LESModel>("LESProperties");

    // Field data
    const label patchI = patch().index();

    const scalarField& muw = lesModel.mu().boundaryField()[patchI];
    const scalarField muSgsw = lesModel.muSgs()().boundaryField()[patchI];

    const scalarField& alphaw = lesModel.alpha().boundaryField()[patchI];
    scalarField& alphaSgsw = *this;

    const fvPatchVectorField& Uw = lesModel.U().boundaryField()[patchI];
    const scalarField magUp = mag(Uw.patchInternalField() - Uw);
    const scalarField magGradUw = mag(Uw.snGrad());

    const scalarField& rhow = lesModel.rho().boundaryField()[patchI];
    const fvPatchScalarField& hw =
        patch().lookupPatchField<volScalarField, scalar>(hsName_);

    const scalarField& ry = patch().deltaCoeffs();

    // Heat flux [W/m2] - lagging alphaSgsw
    const scalarField qDot = (alphaw + alphaSgsw)*hw.snGrad();

    // Populate boundary values
    forAll(alphaSgsw, faceI)
    {
        // Calculate uTau using Newton-Raphson iteration
        scalar uTau =
            sqrt((muSgsw[faceI] + muw[faceI])/rhow[faceI]*magGradUw[faceI]);

        if (uTau > ROOTVSMALL)
        {
            label iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUp[faceI]/uTau, maxExp_);
                scalar fkUu = exp(kUu) - 1.0 - kUu*(1.0 + 0.5*kUu);

                scalar f =
                    - uTau/(ry[faceI]*muw[faceI]/rhow[faceI])
                    + magUp[faceI]/uTau
                    + 1.0/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    - 1.0/(ry[faceI]*muw[faceI]/rhow[faceI])
                    - magUp[faceI]/sqr(uTau)
                    - 1.0/E_*kUu*fkUu/uTau;

                scalar uTauNew = uTau - f/df;
                err = mag((uTau - uTauNew)/uTau);
                uTau = uTauNew;

            } while (uTau>VSMALL && err>tolerance_ && ++iter<maxIters_);

            scalar yPlus = uTau/ry[faceI]/(muw[faceI]/rhow[faceI]);

            // Molecular Prandtl number
            scalar Pr = muw[faceI]/alphaw[faceI];

            // Molecular-to-turbulenbt Prandtl number ratio
            scalar Prat = Pr/Prt_;

            // Thermal sublayer thickness
            scalar P = Psmooth(Prat);
            scalar yPlusTherm = this->yPlusTherm(P, Prat);

            // Evaluate new effective thermal diffusivity
            scalar alphaEff = 0.0;
            if (yPlus < yPlusTherm)
            {
                scalar A = qDot[faceI]*rhow[faceI]*uTau/ry[faceI];
                scalar B = qDot[faceI]*Pr*yPlus;
                scalar C = Pr*0.5*rhow[faceI]*uTau*sqr(magUp[faceI]);
                alphaEff = A/(B + C + VSMALL);
            }
            else
            {
                scalar A = qDot[faceI]*rhow[faceI]*uTau/ry[faceI];
                scalar B = qDot[faceI]*Prt_*(1.0/kappa_*log(E_*yPlus) + P);
                scalar magUc = uTau/kappa_*log(E_*yPlusTherm) - mag(Uw[faceI]);
                scalar C =
                    0.5*rhow[faceI]*uTau
                   *(Prt_*sqr(magUp[faceI]) + (Pr - Prt_)*sqr(magUc));
                alphaEff = A/(B + C + VSMALL);
            }

            // Update turbulent thermal diffusivity
            alphaSgsw[faceI] = max(0.0, alphaEff - alphaw[faceI]);

            if (debug)
            {
                Info<< "    uTau           = " << uTau << nl
                    << "    Pr             = " << Pr << nl
                    << "    Prt            = " << Prt_ << nl
                    << "    qDot           = " << qDot[faceI] << nl
                    << "    yPlus          = " << yPlus << nl
                    << "    yPlusTherm     = " << yPlusTherm << nl
                    << "    alphaEff       = " << alphaEff << nl
                    << "    alphaw         = " << alphaw[faceI] << nl
                    << "    alphaSgsw      = " << alphaSgsw[faceI] << nl
                    << endl;
            }
        }
        else
        {
            alphaSgsw[faceI] = 0.0;
        }
    }
}


void alphaSgsJayatillekeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("Prt") << Prt_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("hs") << hsName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphaSgsJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
