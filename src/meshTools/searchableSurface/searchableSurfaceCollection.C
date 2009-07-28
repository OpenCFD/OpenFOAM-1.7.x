/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "searchableSurfaceCollection.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"
#include "Time.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableSurfaceCollection, 0);
addToRunTimeSelectionTable
(
    searchableSurface,
    searchableSurfaceCollection,
    dict
);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceCollection::findNearest
(
    const pointField& samples,
    scalarField& minDistSqr,
    List<pointIndexHit>& nearestInfo,
    labelList& nearestSurf
) const
{
    // Initialise
    nearestInfo.setSize(samples.size());
    nearestInfo = pointIndexHit();
    nearestSurf.setSize(samples.size());
    nearestSurf = -1;

    List<pointIndexHit> hitInfo(samples.size());

    const scalarField localMinDistSqr(samples.size(), GREAT);

    forAll(subGeom_, surfI)
    {
        subGeom_[surfI].findNearest
        (
            cmptDivide  // Transform then divide
            (
                transform_[surfI].localPosition(samples),
                scale_[surfI]
            ),
            localMinDistSqr,
            hitInfo
        );

        forAll(hitInfo, pointI)
        {
            if (hitInfo[pointI].hit())
            {
                // Rework back into global coordinate sys. Multiply then
                // transform
                point globalPt = transform_[surfI].globalPosition
                (
                    cmptMultiply
                    (
                        hitInfo[pointI].rawPoint(),
                        scale_[surfI]
                    )
                );

                scalar distSqr = magSqr(globalPt - samples[pointI]);

                if (distSqr < minDistSqr[pointI])
                {
                    minDistSqr[pointI] = distSqr;
                    nearestInfo[pointI].setPoint(globalPt);
                    nearestInfo[pointI].setHit();
                    nearestInfo[pointI].setIndex(hitInfo[pointI].index());
                    nearestSurf[pointI] = surfI;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceCollection::searchableSurfaceCollection
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    instance_(dict.size()),
    scale_(dict.size()),
    transform_(dict.size()),
    subGeom_(dict.size()),
    mergeSubRegions_(dict.lookup("mergeSubRegions"))
{
    Info<< "SearchableCollection : " << name() << endl;

    label surfI = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (dict.isDict(iter().keyword()))
        {
            instance_[surfI] = iter().keyword();

            const dictionary& subDict = dict.subDict(instance_[surfI]);

            scale_[surfI] = subDict.lookup("scale");
            transform_.set
            (
                surfI,
                coordinateSystem::New
                (
                    "",
                    subDict.subDict("transform")
                )
            );

            const word subGeomName(subDict.lookup("surface"));
            //Pout<< "Trying to find " << subGeomName << endl;

            const searchableSurface& s =
                io.db().lookupObject<searchableSurface>(subGeomName);

            subGeom_.set(surfI, &const_cast<searchableSurface&>(s));

            Info<< "    instance : " << instance_[surfI] << endl;
            Info<< "    surface  : " << s.name() << endl;
            Info<< "    scale    : " << scale_[surfI] << endl;
            Info<< "    coordsys : " << transform_[surfI] << endl;

            surfI++;
        }
    }
    instance_.setSize(surfI);
    scale_.setSize(surfI);
    transform_.setSize(surfI);
    subGeom_.setSize(surfI);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceCollection::~searchableSurfaceCollection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableSurfaceCollection::regions() const
{
    if (regions_.size() == 0)
    {
        regionOffset_.setSize(subGeom_.size());

        DynamicList<word> allRegions;
        forAll(subGeom_, surfI)
        {
            regionOffset_[surfI] = allRegions.size();

            if (mergeSubRegions_)
            {
                // Single name regardless how many regions subsurface has
                allRegions.append(instance_[surfI] + "_" + Foam::name(surfI));
            }
            else
            {
                const wordList& subRegions = subGeom_[surfI].regions();

                forAll(subRegions, i)
                {
                    allRegions.append(instance_[surfI] + "_" + subRegions[i]);
                }
            }
        }
        regions_.transfer(allRegions.shrink());
    }
    return regions_;
}


Foam::label Foam::searchableSurfaceCollection::size() const
{
    label n = 0;
    forAll(subGeom_, surfI)
    {
        n += subGeom_[surfI].size();
    }
    return n;
}


void Foam::searchableSurfaceCollection::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& nearestInfo
) const
{
    // How to scale distance?
    scalarField minDistSqr(nearestDistSqr);

    labelList nearestSurf;
    findNearest
    (
        samples,
        minDistSqr,
        nearestInfo,
        nearestSurf
    );
}


void Foam::searchableSurfaceCollection::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());
    info = pointIndexHit();

    // Current nearest (to start) intersection
    pointField nearest(end);

    List<pointIndexHit> hitInfo(start.size());

    forAll(subGeom_, surfI)
    {
        // Starting point
        tmp<pointField> e0 = cmptDivide
        (
            transform_[surfI].localPosition
            (
                start
            ),
            scale_[surfI]
        );

        // Current best end point
        tmp<pointField> e1 = cmptDivide
        (
            transform_[surfI].localPosition
            (
                nearest
            ),
            scale_[surfI]
        );

        subGeom_[surfI].findLine(e0, e1, hitInfo);

        forAll(hitInfo, pointI)
        {
            if (hitInfo[pointI].hit())
            {
                // Transform back to global coordinate sys.
                nearest[pointI] = transform_[surfI].globalPosition
                (
                    cmptMultiply
                    (
                        hitInfo[pointI].rawPoint(),
                        scale_[surfI]
                    )
                );
                info[pointI] = hitInfo[pointI];
                info[pointI].rawPoint() = nearest[pointI];
            }
        }
    }


    // Debug check
    if (false)
    {
        forAll(info, pointI)
        {
            if (info[pointI].hit())
            {
                vector n(end[pointI] - start[pointI]);
                scalar magN = mag(n);

                if (magN > SMALL)
                {
                    n /= mag(n);

                    scalar s = ((info[pointI].rawPoint()-start[pointI])&n);

                    if (s < 0 || s > 1)
                    {
                        FatalErrorIn
                        (
                            "searchableSurfaceCollection::findLine(..)"
                        )   << "point:" << info[pointI]
                            << " s:" << s
                            << " outside vector "
                            << " start:" << start[pointI]
                            << " end:" << end[pointI]
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}


void Foam::searchableSurfaceCollection::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    // To be done ...
    findLine(start, end, info);
}


void Foam::searchableSurfaceCollection::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit> >& info
) const
{
    // To be done. Assume for now only one intersection.
    List<pointIndexHit> nearestInfo;
    findLine(start, end, nearestInfo);

    info.setSize(start.size());
    forAll(info, pointI)
    {
        if (nearestInfo[pointI].hit())
        {
            info[pointI].setSize(1);
            info[pointI][0] = nearestInfo[pointI];
        }
        else
        {
            info[pointI].clear();
        }
    }
}


void Foam::searchableSurfaceCollection::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    if (subGeom_.size() == 0)
    {}
    else if (subGeom_.size() == 1)
    {
        if (mergeSubRegions_)
        {
            region.setSize(info.size());
            region = regionOffset_[0];
        }
        else
        {
            subGeom_[0].getRegion(info, region);
        }
    }
    else
    {
        region.setSize(info.size());
        region = -1;

        // Which region did point come from. Retest for now to see which
        // surface it originates from - crap solution! Should use global indices
        // in index inside pointIndexHit to do this better.

        pointField samples(info.size());
        forAll(info, pointI)
        {
            if (info[pointI].hit())
            {
                samples[pointI] = info[pointI].hitPoint();
            }
            else
            {
                samples[pointI] = vector::zero;
            }
        }
        //scalarField minDistSqr(info.size(), SMALL);
        scalarField minDistSqr(info.size(), GREAT);

        labelList nearestSurf;
        List<pointIndexHit> nearestInfo;
        findNearest
        (
            samples,
            minDistSqr,
            nearestInfo,
            nearestSurf
        );

        // Check
        {
            forAll(info, pointI)
            {
                if (info[pointI].hit() && nearestSurf[pointI] == -1)
                {
                    FatalErrorIn
                    (
                        "searchableSurfaceCollection::getRegion(..)"
                    )   << "pointI:" << pointI
                        << " sample:" << samples[pointI]
                        << " nearest:" << nearestInfo[pointI]
                        << " nearestsurf:" << nearestSurf[pointI]
                        << abort(FatalError);
                }
            }
        }

        forAll(subGeom_, surfI)
        {
            // Collect points from my surface
            labelList indices(findIndices(nearestSurf, surfI));

            if (mergeSubRegions_)
            {
                forAll(indices, i)
                {
                    region[indices[i]] = regionOffset_[surfI];
                }
            }
            else
            {
                labelList surfRegion;
                subGeom_[surfI].getRegion
                (
                    UIndirectList<pointIndexHit>(info, indices),
                    surfRegion
                );
                forAll(indices, i)
                {
                    region[indices[i]] = regionOffset_[surfI] + surfRegion[i];
                }
            }
        }
    }
}


void Foam::searchableSurfaceCollection::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    if (subGeom_.size() == 0)
    {}
    else if (subGeom_.size() == 1)
    {
        subGeom_[0].getNormal(info, normal);
    }
    else
    {
        normal.setSize(info.size());

        // See above - crap retest to find surface point originates from.
        pointField samples(info.size());
        forAll(info, pointI)
        {
            if (info[pointI].hit())
            {
                samples[pointI] = info[pointI].hitPoint();
            }
            else
            {
                samples[pointI] = vector::zero;
            }
        }
        //scalarField minDistSqr(info.size(), SMALL);
        scalarField minDistSqr(info.size(), GREAT);

        labelList nearestSurf;
        List<pointIndexHit> nearestInfo;
        findNearest
        (
            samples,
            minDistSqr,
            nearestInfo,
            nearestSurf
        );


        forAll(subGeom_, surfI)
        {
            // Collect points from my surface
            labelList indices(findIndices(nearestSurf, surfI));

            vectorField surfNormal;
            subGeom_[surfI].getNormal
            (
                UIndirectList<pointIndexHit>(info, indices),
                surfNormal
            );
            forAll(indices, i)
            {
                normal[indices[i]] = surfNormal[i];
            }
        }
    }
}


void Foam::searchableSurfaceCollection::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorIn
    (
        "searchableSurfaceCollection::getVolumeType(const pointField&"
        ", List<volumeType>&) const"
    )   << "Volume type not supported for collection."
        << exit(FatalError);
}


// ************************************************************************* //
