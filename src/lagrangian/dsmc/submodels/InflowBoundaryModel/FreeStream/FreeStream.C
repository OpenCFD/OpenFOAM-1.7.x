/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "FreeStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::FreeStream<CloudType>::FreeStream
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InflowBoundaryModel<CloudType>(dict, cloud, typeName),
    patches_(),
    moleculeTypeIds_(),
    numberDensities_(),
    particleFluxAccumulators_()
{
    // Identify which patches to use

    DynamicList<label> patches;

    forAll(cloud.mesh().boundaryMesh(), p)
    {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[p];

        if (patch.type() == polyPatch::typeName)
        {
            patches.append(p);
        }
    }

    patches_.transfer(patches);

    const dictionary& numberDensitiesDict
    (
        this->coeffDict().subDict("numberDensities")
    );

    List<word> molecules(numberDensitiesDict.toc());

    // Initialise the particleFluxAccumulators_
    particleFluxAccumulators_.setSize(patches_.size());

    forAll(patches_, p)
    {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patches_[p]];

        particleFluxAccumulators_[p] = List<Field<scalar> >
        (
            molecules.size(),
            Field<scalar>(patch.size(), 0.0)
        );
    }

    moleculeTypeIds_.setSize(molecules.size());

    numberDensities_.setSize(molecules.size());

    forAll(molecules, i)
    {
        numberDensities_[i] = readScalar
        (
            numberDensitiesDict.lookup(molecules[i])
        );

        moleculeTypeIds_[i] = findIndex(cloud.typeIdList(), molecules[i]);

        if (moleculeTypeIds_[i] == -1)
        {
            FatalErrorIn
            (
                "Foam::FreeStream<CloudType>::FreeStream"
                "("
                    "const dictionary&, "
                    "CloudType&"
                ")"
            )   << "typeId " << molecules[i] << "not defined in cloud." << nl
                << abort(FatalError);
        }
    }

    numberDensities_ /= cloud.nParticle();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::FreeStream<CloudType>::~FreeStream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::FreeStream<CloudType>::inflow()
{
    CloudType& cloud(this->owner());

    const polyMesh& mesh(cloud.mesh());

    const scalar deltaT = cloud.cachedDeltaT();

    Random& rndGen(cloud.rndGen());

    scalar sqrtPi = sqrt(mathematicalConstant::pi);

    label particlesInserted = 0;

    const volScalarField::GeometricBoundaryField& boundaryT
    (
        cloud.T().boundaryField()
    );

    const volVectorField::GeometricBoundaryField& boundaryU
    (
        cloud.U().boundaryField()
    );

    forAll(patches_, p)
    {
        label patchI = patches_[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        // Add mass to the accumulators.  negative face area dotted with the
        // velocity to point flux into the domain.

        // Take a reference to the particleFluxAccumulator for this patch
        List<Field<scalar> >& pFA = particleFluxAccumulators_[p];

        forAll(pFA, i)
        {
            label typeId = moleculeTypeIds_[i];

            scalar mass = cloud.constProps(typeId).mass();

            if (min(boundaryT[patchI]) < SMALL)
            {
                FatalErrorIn ("Foam::FreeStream<CloudType>::inflow()")
                    << "Zero boundary temperature detected, check boundaryT condition." << nl
                    << nl << abort(FatalError);
            }

            scalarField mostProbableSpeed
            (
                cloud.maxwellianMostProbableSpeed
                (
                    boundaryT[patchI],
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal (which points
            // out of the domain, so it must be negated), dividing by the most
            // probable speed to form molecularSpeedRatio * cosTheta

            scalarField sCosTheta =
                boundaryU[patchI]
              & -patch.faceAreas()/mag(patch.faceAreas())
               /mostProbableSpeed;

            // From Bird eqn 4.22

            pFA[i] +=
                mag(patch.faceAreas())*numberDensities_[i]*deltaT
               *mostProbableSpeed
               *(
                   exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
               /(2.0*sqrtPi);
        }

        forAll(patch, f)
        {
            // Loop over all faces as the outer loop to avoid calculating
            // geometrical properties multiple times.

            labelList faceVertices = patch[f];

            label nVertices = faceVertices.size();

            label globalFaceIndex = f + patch.start();

            label cell = mesh.faceOwner()[globalFaceIndex];

            const vector& fC = patch.faceCentres()[f];

            scalar fA = mag(patch.faceAreas()[f]);

            // Cumulative triangle area fractions
            List<scalar> cTriAFracs(nVertices);

            for (label v = 0; v < nVertices - 1; v++)
            {
                const point& vA = mesh.points()[faceVertices[v]];

                const point& vB = mesh.points()[faceVertices[(v + 1)]];

                cTriAFracs[v] =
                    0.5*mag((vA - fC)^(vB - fC))/fA
                  + cTriAFracs[max((v - 1), 0)];
            }

            // Force the last area fraction value to 1.0 to avoid any
            // rounding/non-flat face errors giving a value < 1.0
            cTriAFracs[nVertices - 1] = 1.0;

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = patch.faceAreas()[f];
            n /= -mag(n);

            // Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - (mesh.points()[faceVertices[0]]);
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            scalar faceTemperature = boundaryT[patchI][f];

            const vector& faceVelocity = boundaryU[patchI][f];

            forAll(pFA, i)
            {
                scalar& faceAccumulator = pFA[i][f];

                // Number of whole particles to insert
                label nI = max(label(faceAccumulator), 0);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of faceAccumulator
                if ((faceAccumulator - nI) > rndGen.scalar01())
                {
                    nI++;
                }

                faceAccumulator -= nI;

                label typeId = moleculeTypeIds_[i];

                scalar mass = cloud.constProps(typeId).mass();

                for (label i = 0; i < nI; i++)
                {
                    // Choose a triangle to insert on, based on their relative
                    // area

                    scalar triSelection = rndGen.scalar01();

                    // Selected triangle
                    label sTri = -1;

                    forAll(cTriAFracs, tri)
                    {
                        sTri = tri;

                        if (cTriAFracs[tri] >= triSelection)
                        {
                            break;
                        }
                    }

                    // Randomly distribute the points on the triangle, using the
                    // method from:
                    // Generating Random Points in Triangles
                    // by Greg Turk
                    // from "Graphics Gems", Academic Press, 1990
                    // http://tog.acm.org/GraphicsGems/gems/TriPoints.c

                    const point& A = fC;
                    const point& B = mesh.points()[faceVertices[sTri]];
                    const point& C =
                    mesh.points()[faceVertices[(sTri + 1) % nVertices]];

                    scalar s = rndGen.scalar01();
                    scalar t = sqrt(rndGen.scalar01());

                    point p = (1 - t)*A + (1 - s)*t*B + s*t*C;

                    // Velocity generation

                    scalar mostProbableSpeed
                    (
                        cloud.maxwellianMostProbableSpeed
                        (
                            faceTemperature,
                            mass
                        )
                    );

                    scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                    // Coefficients required for Bird eqn 12.5
                    scalar uNormProbCoeffA =
                        sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                    scalar uNormProbCoeffB =
                        0.5*
                        (
                            1.0
                          + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                        );

                    // Equivalent to the QA value in Bird's DSMC3.FOR
                    scalar randomScaling = 3.0;

                    if (sCosTheta < -3)
                    {
                        randomScaling = mag(sCosTheta) + 1;
                    }

                    scalar P = -1;

                    // Normalised candidates for the normal direction velocity
                    // component
                    scalar uNormal;
                    scalar uNormalThermal;

                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.scalar01() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.scalar01());

                    vector U =
                        sqrt(CloudType::kb*faceTemperature/mass)
                       *(
                            rndGen.GaussNormal()*t1
                          + rndGen.GaussNormal()*t2
                        )
                      + (t1 & faceVelocity)*t1
                      + (t2 & faceVelocity)*t2
                      + mostProbableSpeed*uNormal*n;

                    scalar Ei = cloud.equipartitionInternalEnergy
                    (
                        faceTemperature,
                        cloud.constProps(typeId).internalDegreesOfFreedom()
                    );

                    cloud.addNewParcel
                    (
                        p,
                        U,
                        Ei,
                        cell,
                        typeId
                    );

                    particlesInserted++;
                }
            }
        }
    }

    reduce(particlesInserted, sumOp<label>());

    Info<< "    Particles inserted              = "
        << particlesInserted << endl;

}


// ************************************************************************* //
