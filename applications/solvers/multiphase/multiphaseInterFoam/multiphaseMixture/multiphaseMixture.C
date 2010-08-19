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

#include "multiphaseMixture.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "Time.H"
#include "subCycle.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const scalar Foam::multiphaseMixture::convertToRad =
    Foam::mathematicalConstant::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseMixture::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        alphas_ += level*iter();
        level += 1.0;
    }

    alphas_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseMixture::multiphaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    transportModel(U, phi),
    phases_(lookup("phases"), phase::iNew(U, phi)),
    refPhase_(*phases_.lookup(word(lookup("refPhase")))),

    mesh_(U.mesh()),
    U_(U),
    phi_(phi),

    rhoPhi_
    (
        IOobject
        (
            "rho*phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("rho*phi", dimMass/dimTime, 0.0)
    ),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphas", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
{
    calcAlphas();
    alphas_.write();

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        alphaTable_.add(iter());
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::rho() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> trho = iter()*iter().rho();

    for(++iter; iter != phases_.end(); ++iter)
    {
        trho() += iter()*iter().rho();
    }

    return trho;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::mu() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tmu = iter()*iter().rho()*iter().nu();

    for(++iter; iter != phases_.end(); ++iter)
    {
        tmu() += iter()*iter().rho()*iter().nu();
    }

    return tmu;
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseMixture::muf() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<surfaceScalarField> tmuf =
        fvc::interpolate(iter())*iter().rho()*fvc::interpolate(iter().nu());

    for(++iter; iter != phases_.end(); ++iter)
    {
        tmuf() +=
            fvc::interpolate(iter())*iter().rho()*fvc::interpolate(iter().nu());
    }

    return tmuf;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseMixture::nuf() const
{
    return muf()/fvc::interpolate(rho());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::surfaceTensionForce() const
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "surfaceTensionForce",
                dimensionSet(1, -2, -2, 0, 0),
                0.0
            )
        )
    );

    surfaceScalarField& stf = tstf();

    forAllConstIter(PtrDictionary<phase>, phases_, iter1)
    {
        const phase& alpha1 = iter1();

        PtrDictionary<phase>::const_iterator iter2 = iter1;
        ++iter2;

        for(; iter2 != phases_.end(); ++iter2)
        {
            const phase& alpha2 = iter2();

            sigmaTable::const_iterator sigma =
                sigmas_.find(interfacePair(alpha1, alpha2));

            if (sigma == sigmas_.end())
            {
                FatalErrorIn("multiphaseMixture::surfaceTensionForce() const")
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar("sigma", dimSigma_, sigma())
               *fvc::interpolate(K(alpha1, alpha2))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


void Foam::multiphaseMixture::solve()
{
    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        iter().correct();
    }

    const Time& runTime = mesh_.time();

    label nAlphaSubCycles
    (
        readLabel
        (
            mesh_.solutionDict().subDict("PISO").lookup("nAlphaSubCycles")
        )
    );

    label nAlphaCorr
    (
        readLabel(mesh_.solutionDict().subDict("PISO").lookup("nAlphaCorr"))
    );

    bool cycleAlpha
    (
        Switch(mesh_.solutionDict().subDict("PISO").lookup("cycleAlpha"))
    );

    scalar cAlpha
    (
        readScalar(mesh_.solutionDict().subDict("PISO").lookup("cAlpha"))
    );


    volScalarField& alpha = phases_.first();

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum = 0.0*rhoPhi_;
        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas(nAlphaCorr, cycleAlpha, cAlpha);
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphas(nAlphaCorr, cycleAlpha, cAlpha);
    }
}


void Foam::multiphaseMixture::correct()
{}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseMixture::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf =
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2));

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseMixture::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::multiphaseMixture::correctContactAngle
(
    const phase& alpha1,
    const phase& alpha2,
    surfaceVectorField::GeometricBoundaryField& nHatb
) const
{
    const volScalarField::GeometricBoundaryField& gbf
        = refPhase_.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& acap =
                refCast<const alphaContactAngleFvPatchScalarField>(gbf[patchi]);

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch =
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi];

            alphaContactAngleFvPatchScalarField::thetaPropsTable::
                const_iterator tp =
                acap.thetaProps().find(interfacePair(alpha1, alpha2));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorIn
                (
                    "multiphaseMixture::correctContactAngle"
                    "(const phase& alpha1, const phase& alpha2, "
                    "fvPatchVectorFieldField& nHatb) const"
                )   << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            bool matched = (tp.key().first() == alpha1.name());

            scalar theta0 = convertToRad*tp().theta0(matched);
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > SMALL)
            {
                scalar thetaA = convertToRad*tp().thetaA(matched);
                scalar thetaR = convertToRad*tp().thetaR(matched);

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall =
                    U_.boundaryField()[patchi].patchInternalField()
                  - U_.boundaryField()[patchi];
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall =
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch;

                // Normalise nWall
                nWall /= (mag(nWall) + SMALL);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall = nWall & Uwall;

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            scalarField a12 = nHatPatch & AfHatPatch;

            scalarField b1 = cos(theta);

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det = 1.0 - a12*a12;

            scalarField a = (b1 - a12*b2)/det;
            scalarField b = (b2 - a12*b1)/det;

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::K
(
    const phase& alpha1,
    const phase& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle(alpha1, alpha2, tnHatfv().boundaryField());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::nearInterface() const
{
    tmp<surfaceScalarField> tnearInt
    (
        new surfaceScalarField
        (
            IOobject
            (
                "nearInterface",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("nearInterface", dimless, 0.0)
        )
    );

    forAllConstIter(PtrDictionary<phase>, phases_, iter)
    {
        surfaceScalarField alphaf = fvc::interpolate(iter());
        tnearInt() = max(tnearInt(), pos(alphaf - 0.01)*pos(0.99 - alphaf));
    }

    return tnearInt;
}


void Foam::multiphaseMixture::solveAlphas
(
    const label nAlphaCorr,
    const bool cycleAlpha,
    const scalar cAlpha
)
{
    static label nSolves=-1;
    nSolves++;

    word alphaScheme("div(phi,alpha)");
    word alphacScheme("div(phirb,alpha)");

    tmp<fv::convectionScheme<scalar> > mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh_,
            alphaTable_,
            phi_,
            mesh_.divScheme(alphaScheme)
        )
    );

    surfaceScalarField phic = mag(phi_/mesh_.magSf());
    phic = min(cAlpha*phic, max(phic));

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        phase* refPhasePtr = &refPhase_;

        if (cycleAlpha)
        {
            PtrDictionary<phase>::iterator refPhaseIter = phases_.begin();
            for(label i=0; i<nSolves%phases_.size(); i++)
            {
                ++refPhaseIter;
            }
            refPhasePtr = &refPhaseIter();
        }

        phase& refPhase = *refPhasePtr;

        volScalarField refPhaseNew = refPhase;
        refPhaseNew == 1.0;

        rhoPhi_ = phi_*refPhase.rho();

        forAllIter(PtrDictionary<phase>, phases_, iter)
        {
            phase& alpha = iter();

            if (&alpha == &refPhase) continue;

            fvScalarMatrix alphaEqn
            (
                fvm::ddt(alpha)
              + mvConvection->fvmDiv(phi_, alpha)
            );

            forAllIter(PtrDictionary<phase>, phases_, iter2)
            {
                phase& alpha2 = iter2();

                if (&alpha2 == &alpha) continue;

                surfaceScalarField phir = phic*nHatf(alpha, alpha2);
                surfaceScalarField phirb12 =
                    -fvc::flux(-phir, alpha2, alphacScheme);

                alphaEqn += fvm::div(phirb12, alpha, alphacScheme);
            }

            alphaEqn.solve(mesh_.solver("alpha"));

            rhoPhi_ += alphaEqn.flux()*(alpha.rho() - refPhase.rho());

            Info<< alpha.name() << " volume fraction, min, max = "
                << alpha.weightedAverage(mesh_.V()).value()
                << ' ' << min(alpha).value()
                << ' ' << max(alpha).value()
                << endl;

            refPhaseNew == refPhaseNew - alpha;
        }

        refPhase == refPhaseNew;

        Info<< refPhase.name() << " volume fraction, min, max = "
            << refPhase.weightedAverage(mesh_.V()).value()
            << ' ' << min(refPhase).value()
            << ' ' << max(refPhase).value()
            << endl;
    }

    calcAlphas();
}


bool Foam::multiphaseMixture::read()
{
    if (transportModel::read())
    {
        bool readOK = true;

        PtrList<entry> phaseData(lookup("phases"));
        label phasei = 0;

        forAllIter(PtrDictionary<phase>, phases_, iter)
        {
            readOK &= iter().read(phaseData[phasei++].dict());
        }

        lookup("sigmas") >> sigmas_;

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
