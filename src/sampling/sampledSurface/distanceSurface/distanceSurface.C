/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "distanceSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "isoSurface.H"
// #include "isoSurfaceCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distanceSurface, 0);
    addToRunTimeSelectionTable(sampledSurface, distanceSurface, word);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distanceSurface::createGeometry()
{
    if (debug)
    {
        Pout<< "distanceSurface::createGeometry :updating geometry." << endl;
    }

    // Clear any stored topologies
    facesPtr_.clear();

    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvm,
            dimensionedScalar("zero", dimLength, 0)
        )
    );
    volScalarField& cellDistance = cellDistancePtr_();

    // Internal field
    {
        const pointField& cc = fvm.C();
        scalarField& fld = cellDistance.internalField();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            cc,
            scalarField(cc.size(), GREAT),
            nearest
        );

        if (signed_)
        {
            vectorField normal;
            surfPtr_().getNormal(nearest, normal);

            forAll(nearest, i)
            {
                vector d(cc[i]-nearest[i].hitPoint());

                if ((d&normal[i]) > 0)
                {
                    fld[i] = Foam::mag(d);
                }
                else
                {
                    fld[i] = -Foam::mag(d);
                }
            }
        }
        else
        {
            forAll(nearest, i)
            {
                fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
            }
        }
    }

    // Patch fields
    {
        forAll(fvm.C().boundaryField(), patchI)
        {
            const pointField& cc = fvm.C().boundaryField()[patchI];
            fvPatchScalarField& fld = cellDistance.boundaryField()[patchI];

            List<pointIndexHit> nearest;
            surfPtr_().findNearest
            (
                cc,
                scalarField(cc.size(), GREAT),
                nearest
            );

            if (signed_)
            {
                vectorField normal;
                surfPtr_().getNormal(nearest, normal);

                forAll(nearest, i)
                {
                    vector d(cc[i]-nearest[i].hitPoint());

                    if ((d&normal[i]) > 0)
                    {
                        fld[i] = Foam::mag(d);
                    }
                    else
                    {
                        fld[i] = -Foam::mag(d);
                    }
                }
            }
            else
            {
                forAll(nearest, i)
                {
                    fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
                }
            }
        }
    }


    // On processor patches the mesh.C() will already be the cell centre
    // on the opposite side so no need to swap cellDistance.


    // Distance to points
    pointDistance_.setSize(fvm.nPoints());
    {
        const pointField& pts = fvm.points();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            pts,
            scalarField(pts.size(), GREAT),
            nearest
        );

        if (signed_)
        {
            vectorField normal;
            surfPtr_().getNormal(nearest, normal);

            forAll(nearest, i)
            {
                vector d(pts[i]-nearest[i].hitPoint());

                if ((d&normal[i]) > 0)
                {
                    pointDistance_[i] = Foam::mag(d);
                }
                else
                {
                    pointDistance_[i] = -Foam::mag(d);
                }
            }
        }
        else
        {
            forAll(nearest, i)
            {
                pointDistance_[i] = Foam::mag(pts[i]-nearest[i].hitPoint());
            }
        }
    }


    if (debug)
    {
        Pout<< "Writing cell distance:" << cellDistance.objectPath() << endl;
        cellDistance.write();
        pointScalarField pDist
        (
            IOobject
            (
                "pointDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(fvm),
            dimensionedScalar("zero", dimLength, 0)
        );
        pDist.internalField() = pointDistance_;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }


    //- Direct from cell field and point field.
    isoSurfPtr_.reset
    (
        new isoSurface
        (
            cellDistance,
            pointDistance_,
            distance_,
            regularise_
        )
    );

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distanceSurface::distanceSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.lookup("surfaceType"),
            IOobject
            (
                dict.lookupOrDefault("surfaceName", name),  // name
                mesh.time().constant(),                     // directory
                "triSurface",                               // instance
                mesh.time(),                                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    distance_(readScalar(dict.lookup("distance"))),
    signed_(readBool(dict.lookup("signed"))),
    regularise_(dict.lookupOrDefault("regularise", true)),
    zoneName_(word::null),
    needsUpdate_(true),
    isoSurfPtr_(NULL),
    facesPtr_(NULL)
{
//    dict.readIfPresent("zone", zoneName_);
//
//    if (debug && zoneName_.size())
//    {
//        if (mesh.cellZones().findZoneID(zoneName_) < 0)
//        {
//            Info<< "cellZone \"" << zoneName_
//                << "\" not found - using entire mesh" << endl;
//        }
//    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distanceSurface::~distanceSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::distanceSurface::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::distanceSurface::expire()
{
    if (debug)
    {
        Pout<< "distanceSurface::expire :"
            << " have-facesPtr_:" << facesPtr_.valid()
            << " needsUpdate_:" << needsUpdate_ << endl;
    }

    // Clear any stored topologies
    facesPtr_.clear();

    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::distanceSurface::update()
{
    if (debug)
    {
        Pout<< "distanceSurface::update :"
            << " have-facesPtr_:" << facesPtr_.valid()
            << " needsUpdate_:" << needsUpdate_ << endl;
    }

    if (!needsUpdate_)
    {
        return false;
    }

    createGeometry();

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField>
Foam::distanceSurface::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::distanceSurface::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::distanceSurface::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::distanceSurface::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::distanceSurface::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::distanceSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::distanceSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::distanceSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::distanceSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::distanceSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::distanceSurface::print(Ostream& os) const
{
    os  << "distanceSurface: " << name() << " :"
        << "  surface:" << surfPtr_().name()
        << "  distance:" << distance_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
