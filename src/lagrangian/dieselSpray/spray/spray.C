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

#include "spray.H"

#include "atomizationModel.H"
#include "breakupModel.H"
#include "collisionModel.H"
#include "dispersionModel.H"
#include "dragModel.H"
#include "evaporationModel.H"
#include "heatTransferModel.H"
#include "injectorModel.H"
#include "wallModel.H"

#include "basicMultiComponentMixture.H"

#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(IOPtrList<injector>, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::spray::spray
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& p,
    const volScalarField& T,
    const basicMultiComponentMixture& composition,
    const PtrList<gasThermoPhysics>& gasProperties,
    const dictionary&,
    const dimensionedVector& g,
    bool readFields
)
:
    Cloud<parcel>(U.mesh(), false), // suppress className checking on positions
    runTime_(U.time()),
    time0_(runTime_.value()),
    mesh_(U.mesh()),
    rndGen_(label(0)),
    g_(g.value()),

    U_(U),
    rho_(rho),
    p_(p),
    T_(T),

    sprayProperties_
    (
        IOobject
        (
            "sprayProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    ambientPressure_(p_.average().value()),
    ambientTemperature_(T_.average().value()),

    injectors_
    (
        IOobject
        (
            "injectorProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        injector::iNew(U.time())
    ),
    atomization_
    (
        atomizationModel::New
        (
            sprayProperties_,
            *this
        )
    ),
    drag_
    (
        dragModel::New
        (
            sprayProperties_
        )
    ),
    evaporation_
    (
        evaporationModel::New
        (
            sprayProperties_
        )
    ),
    heatTransfer_
    (
        heatTransferModel::New
        (
            sprayProperties_
        )
    ),
    wall_
    (
        wallModel::New
        (
            sprayProperties_,
            U,
            *this
        )
    ),
    breakupModel_
    (
        breakupModel::New
        (
            sprayProperties_,
            *this
        )
    ),
    collisionModel_
    (
        collisionModel::New
        (
            sprayProperties_,
            *this,
            rndGen_
        )
    ),
    dispersionModel_
    (
        dispersionModel::New
        (
            sprayProperties_,
            *this
        )
    ),

    fuels_
    (
        liquidMixture::New
        (
            mesh_.lookupObject<dictionary>("thermophysicalProperties")
        )
    ),
    injectorModel_
    (
        injectorModel::New
        (
            sprayProperties_,
            *this
        )
    ),

    subCycles_(readLabel(sprayProperties_.lookup("subCycles"))),

    gasProperties_(gasProperties),
    composition_(composition),

    liquidToGasIndex_(fuels_->components().size(), -1),
    gasToLiquidIndex_(composition.Y().size(), -1),
    isLiquidFuel_(composition.Y().size(), false),

    twoD_(0),
    axisOfSymmetry_(vector::zero),
    axisOfWedge_(vector(0,0,0)),
    axisOfWedgeNormal_(vector(0,0,0)),
    angleOfWedge_(0.0),

    interpolationSchemes_(sprayProperties_.subDict("interpolationSchemes")),
    UInterpolator_(NULL),
    rhoInterpolator_(NULL),
    pInterpolator_(NULL),
    TInterpolator_(NULL),

    sms_(mesh_.nCells(), vector::zero),
    shs_(mesh_.nCells(), 0.0),
    srhos_(fuels_->components().size()),

    totalInjectedLiquidMass_(0.0),
    injectedLiquidKE_(0.0)

{
    // create the evaporation source fields
    forAll(srhos_, i)
    {
        srhos_.set(i, new scalarField(mesh_.nCells(), 0.0));
    }

    // Write some information about injection parameters
    forAll(injectors_, i)
    {
        const injectorType& it = injectors_[i].properties();

        scalar v = injection().averageVelocity(i);

        scalar ip = it.integrateTable(it.injectionPressureProfile());
        scalar dt = it.teoi() - it.tsoi();
        Info<< "Average Velocity for injector " << i << ": " << v << " m/s"
            << ", injection pressure = "
            << 1.0e-5*ip/dt << " bar"
            << endl;
    }

    // Check if the case is 2D wedge
    const polyBoundaryMesh& bMesh = mesh().boundaryMesh();
    bool symPlaneExist = false;
    bool wedgeExist = false;
    label patches[2];
    label n=0;

    // check for the type of boundary condition
    forAll(bMesh, patchi)
    {
        if (isA<symmetryPolyPatch>(bMesh[patchi]))
        {
            symPlaneExist = true;
        }
        else if (isA<wedgePolyPatch>(bMesh[patchi]))
        {
            wedgeExist = true;
            patches[n++] = patchi;
        }
    }

    // if wedge exist we assume that this is a 2D run.
    twoD_ = wedgeExist;

    if (twoD_)
    {
        if (n<2)
        {
            FatalErrorIn
            (
                "spray::spray(const volVectorField& U, "
                "const volScalarField& rho, const volScalarField& p, "
                "const volScalarField& T, const combustionMixture& composition,"
                "const PtrList<gasThermoPhsyics>& gaseousFuelProperties, "
                "const dictionary& thermophysicalProperties, "
                "const dimensionedScalar& g)"
            )   << "spray::(...) only one wedgePolyPatch found. "
                   "Please check you BC-setup."
                << abort(FatalError);
        }

        Info<< "Constructing two dimensional spray injection.";

        vector v1 = bMesh[patches[0]].faceAreas()[0];
        vector v2 = bMesh[patches[1]].faceAreas()[0];
        v1 /= mag(v1);
        v2 /= mag(v2);
        axisOfSymmetry_ = v1 ^ v2;
        axisOfSymmetry_ /= mag(axisOfSymmetry_);

        // assuming that 'v2' is the 'front' face
        axisOfWedge_ = axisOfSymmetry_ ^ v2;
        axisOfWedge_ /= mag(axisOfWedge_);

        axisOfWedgeNormal_ = axisOfSymmetry_ ^ axisOfWedge_;
        axisOfWedgeNormal_ /= mag(axisOfWedgeNormal_);

        scalar arcCos = (v1 & v2)/mag(v1);
        angleOfWedge_ = mathematicalConstant::pi - acos(arcCos);

        Info<< "Calculated angle of wedge is "
            << angleOfWedge_*180/mathematicalConstant::pi << " deg."
            << endl;
    }
    else
    {
        if (symPlaneExist)
        {
            angleOfWedge_ = mathematicalConstant::pi;
            Info<< "Constructing 180 deg three dimensional spray injection."
                << endl;
        }
        else
        {
            Info<< "Constructing three dimensional spray injection." << endl;
        }

    }

    // find index mapping between liquid indeces and gas indeces
    label Ns = composition_.Y().size();

    forAll(fuels_->components(), i)
    {
        word liquidName(fuels_->components()[i]);

        for (label j=0; j<Ns; j++)
        {
            word specieName(composition_.Y()[j].name());

            if (specieName == liquidName)
            {
                liquidToGasIndex_[i] = j;
                gasToLiquidIndex_[j] = i;
                isLiquidFuel_[j] = true;
            }
        }
        if (liquidToGasIndex_[i] == -1)
        {
            Info << "In composition:" << endl;
            for (label k=0; k<Ns; k++)
            {
                word specieName(composition_.Y()[k].name());
                Info << specieName << endl;
            }

            FatalError<<
                "The liquid component " << liquidName
                << " does not exist in the species composition.Y() list.\n"
                << "(Probably not defined in <chem.inp>)"
                << abort(FatalError);
        }
    }

    if (readFields)
    {
        parcel::readFields(*this);
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::spray::~spray()
{}


// ************************************************************************* //
