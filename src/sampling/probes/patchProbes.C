/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "patchProbes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOmanip.H"
// For 'nearInfo' helper class only
#include "directMappedPatchBase.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchProbes, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchProbes::findElements(const fvMesh& mesh)
{

     // All the info for nearest. Construct to miss
    List<directMappedPatchBase::nearInfo> nearest(probeLocations_.size());

    if (elementList_.empty())
    {

        elementList_.setSize(probeLocations_.size());

        // Octree based search engine
        meshSearch meshSearchEngine(mesh, false);

        forAll(probeLocations_, probeI)
        {
            const point sample = probeLocations_[probeI];

            label faceI = meshSearchEngine.findNearestBoundaryFace(sample);

            if (faceI == -1)
            {
                nearest[probeI].second().first() = Foam::sqr(GREAT);
                nearest[probeI].second().second() = Pstream::myProcNo();
            }
            else
            {
                const point& fc = mesh.faceCentres()[faceI];
                nearest[probeI].first() = pointIndexHit
                (
                    true,
                    fc,
                    faceI
                );
                nearest[probeI].second().first() = magSqr(fc-sample);
                nearest[probeI].second().second() = Pstream::myProcNo();
            }
        }
    }

    // Find nearest.
    Pstream::listCombineGather(nearest, directMappedPatchBase::nearestEqOp());
    Pstream::listCombineScatter(nearest);

    if (debug)
    {
        Info<< "patchProbes::findElements" << " : " << endl;
        forAll(nearest, sampleI)
        {
            label procI = nearest[sampleI].second().second();
            label localI = nearest[sampleI].first().index();

            Info<< "    " << sampleI << " coord:"<< probeLocations_[sampleI]
                << " found on processor:" << procI
                << " in local cell/face:" << localI
                << " with cc:" << nearest[sampleI].first().rawPoint() << endl;
        }
    }



    // Check if all patchProbes have been found.
    forAll(nearest, sampleI)
    {
        label localI = nearest[sampleI].first().index();

        if (localI == -1)
        {
             if (Pstream::master())
             {
                WarningIn("patchProbes::findElements()")
                    << "Did not find location "
                    <<  nearest[sampleI].second().first()
                    << " in any cell. Skipping location." << endl;
             }
        }
        else
        {
            elementList_[sampleI] = localI;
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchProbes::patchProbes
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    probes(name, obr, dict, loadFromFiles)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchProbes::~patchProbes()
{}


void Foam::patchProbes::write()
{
    if (probeLocations_.size() && checkFieldTypes())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);
    }
}

void Foam::patchProbes::read(const dictionary& dict)
{
    probes::read(dict);
}


// ************************************************************************* //
