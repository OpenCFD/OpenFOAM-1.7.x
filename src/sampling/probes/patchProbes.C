/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011 OpenCFD Ltd.
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
#include "IOmanip.H"
// For 'nearInfo' helper class only
#include "directMappedPatchBase.H"
#include "meshSearch.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchProbes, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchProbes::findElements(const fvMesh& mesh)
{

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    label patchI = bm.findPatchID(patchName_);

    if (patchI == -1)
    {
        FatalErrorIn
        (
            " Foam::patchProbes::findElements(const fvMesh&)"
        )   << " Unknown patch name "
            << patchName_ << endl
            << exit(FatalError);
    }


    // All the info for nearest. Construct to miss
    List<directMappedPatchBase::nearInfo> nearest(probeLocations_.size());

    const polyPatch& pp = bm[patchI];

    if (pp.size() > 0)
    {
        labelList bndFaces(pp.size());
        forAll(bndFaces, i)
        {
            bndFaces[i] =  pp.start() + i;
        }

        treeBoundBox overallBb(pp.points());
        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        const indexedOctree<treeDataFace> boundaryTree
        (
            treeDataFace    // all information needed to search faces
            (
                false,                      // do not cache bb
                mesh,
                bndFaces                    // patch faces only
            ),
            overallBb,                      // overall search domain
            8,                              // maxLevel
            10,                             // leafsize
            3.0                             // duplicity
        );


        if (elementList_.empty())
        {
            elementList_.setSize(probeLocations_.size());

            // Octree based search engine
            //meshSearch meshSearchEngine(mesh, false);

            forAll(probeLocations_, probeI)
            {
                const point sample = probeLocations_[probeI];

                scalar span = boundaryTree.bb().mag();

                pointIndexHit info = boundaryTree.findNearest
                (
                    sample,
                    Foam::sqr(span)
                );

                if (!info.hit())
                {
                    info = boundaryTree.findNearest
                    (
                        sample,
                        Foam::sqr(GREAT)
                    );
                }

                label faceI = boundaryTree.shapes().faceLabels()[info.index()];

                const label patchi = bm.whichPatch(faceI);

                if (isA<emptyPolyPatch>(bm[patchi]))
                {
                    WarningIn
                    (
                        " Foam::patchProbes::findElements(const fvMesh&)"
                    )
                    << " The sample point: " << sample
                    << " belongs to " << patchi
                    << " which is an empty patch. This is not permitted. "
                    << " This sample will not be included "
                    << endl;
                }
                else
                {
                    const point& fc = mesh.faceCentres()[faceI];

                    directMappedPatchBase::nearInfo sampleInfo;

                    sampleInfo.first() = pointIndexHit
                    (
                        true,
                        fc,
                        faceI
                    );

                    sampleInfo.second().first() = magSqr(fc-sample);
                    sampleInfo.second().second() = Pstream::myProcNo();

                    nearest[probeI]= sampleInfo;
                }
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
                << " with cc:" << nearest[sampleI].first().rawPoint()
                << " in patch : "<< pp.name() << endl;
        }
    }

    // Check if all patchProbes have been found.
    forAll(probeLocations_, sampleI)
    {
        label localI = -1;
        if (nearest[sampleI].second().second() == Pstream::myProcNo())
        {
            localI = nearest[sampleI].first().index();
        }

        if (elementList_.empty())
        {
             elementList_.setSize(probeLocations_.size());
        }

        elementList_[sampleI] = localI;
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
    dict.lookup("patchName") >> patchName_;
    probes::read(dict);
}


// ************************************************************************* //
