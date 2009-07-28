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

#include "searchableSurface.H"
#include "shellSurfaces.H"
#include "boundBox.H"
#include "triSurfaceMesh.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "orientedSurface.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{

template<>
const char*
NamedEnum<shellSurfaces::refineMode, 3>::
names[] =
{
    "inside",
    "outside",
    "distance"
};

const NamedEnum<shellSurfaces::refineMode, 3> shellSurfaces::refineModeNames_;

} // End namespace Foam



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::shellSurfaces::setAndCheckLevels
(
    const label shellI,
    const List<Tuple2<scalar, label> >& distLevels
)
{
    if (modes_[shellI] != DISTANCE && distLevels.size() != 1)
    {
        FatalErrorIn
        (
            "shellSurfaces::shellSurfaces"
            "(const searchableSurfaces&, const dictionary&)"
        )   << "For refinement mode "
            << refineModeNames_[modes_[shellI]]
            << " specify only one distance+level."
            << " (its distance gets discarded)"
            << exit(FatalError);
    }
    // Extract information into separate distance and level
    distances_[shellI].setSize(distLevels.size());
    levels_[shellI].setSize(distLevels.size());

    forAll(distLevels, j)
    {
        distances_[shellI][j] = distLevels[j].first();
        levels_[shellI][j] = distLevels[j].second();

        // Check in incremental order
        if (j > 0)
        {
            if
            (
                (distances_[shellI][j] <= distances_[shellI][j-1])
             || (levels_[shellI][j] > levels_[shellI][j-1])
            )
            {
                FatalErrorIn
                (
                    "shellSurfaces::shellSurfaces"
                    "(const searchableSurfaces&, const dictionary&)"
                )   << "For refinement mode "
                    << refineModeNames_[modes_[shellI]]
                    << " : Refinement should be specified in order"
                    << " of increasing distance"
                    << " (and decreasing refinement level)." << endl
                    << "Distance:" << distances_[shellI][j]
                    << " refinementLevel:" << levels_[shellI][j]
                    << exit(FatalError);
            }
        }
    }

    const searchableSurface& shell = allGeometry_[shells_[shellI]];

    if (modes_[shellI] == DISTANCE)
    {
        Info<< "Refinement level according to distance to "
            << shell.name() << endl;
        forAll(levels_[shellI], j)
        {
            Info<< "    level " << levels_[shellI][j]
                << " for all cells within " << distances_[shellI][j]
                << " meter." << endl;
        }
    }
    else
    {
        if (!allGeometry_[shells_[shellI]].hasVolumeType())
        {
            FatalErrorIn
            (
                "shellSurfaces::shellSurfaces"
                "(const searchableSurfaces&"
                ", const PtrList<dictionary>&)"
            )   << "Shell " << shell.name()
                << " does not support testing for "
                << refineModeNames_[modes_[shellI]] << endl
                << "Probably it is not closed."
                << exit(FatalError);
        }

        if (modes_[shellI] == INSIDE)
        {
            Info<< "Refinement level " << levels_[shellI][0]
                << " for all cells inside " << shell.name() << endl;
        }
        else
        {
            Info<< "Refinement level " << levels_[shellI][0]
                << " for all cells outside " << shell.name() << endl;
        }
    }
}


// Specifically orient triSurfaces using a calculated point outside.
// Done since quite often triSurfaces not of consistent orientation which
// is (currently) necessary for sideness calculation
void Foam::shellSurfaces::orient()
{
    // Determine outside point.
    boundBox overallBb = boundBox::invertedBox;

    bool hasSurface = false;

    forAll(shells_, shellI)
    {
        const searchableSurface& s = allGeometry_[shells_[shellI]];

        if (modes_[shellI] != DISTANCE && isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& shell = refCast<const triSurfaceMesh>(s);

            if (shell.triSurface::size())
            {
                const pointField& points = shell.points();

                hasSurface = true;

                boundBox shellBb(points[0], points[0]);
                // Assume surface is compact!
                for (label i = 0; i < points.size(); i++)
                {
                    const point& pt = points[i];
                    shellBb.min() = min(shellBb.min(), pt);
                    shellBb.max() = max(shellBb.max(), pt);
                }

                overallBb.min() = min(overallBb.min(), shellBb.min());
                overallBb.max() = max(overallBb.max(), shellBb.max());
            }
        }
    }

    if (hasSurface)
    {
        const point outsidePt = overallBb.max() + overallBb.span();

        //Info<< "Using point " << outsidePt << " to orient shells" << endl;

        forAll(shells_, shellI)
        {
            const searchableSurface& s = allGeometry_[shells_[shellI]];

            if (modes_[shellI] != DISTANCE && isA<triSurfaceMesh>(s))
            {
                triSurfaceMesh& shell = const_cast<triSurfaceMesh&>
                (
                    refCast<const triSurfaceMesh>(s)
                );

                // Flip surface so outsidePt is outside.
                bool anyFlipped = orientedSurface::orient
                (
                    shell,
                    outsidePt,
                    true
                );

                if (anyFlipped)
                {
                    // orientedSurface will have done a clearOut of the surface.
                    // we could do a clearout of the triSurfaceMeshes::trees()
                    // but these aren't affected by orientation
                    // (except for cached
                    // sideness which should not be set at this point.
                    // !!Should check!)

                    Info<< "shellSurfaces : Flipped orientation of surface "
                        << s.name()
                        << " so point " << outsidePt << " is outside." << endl;
                }
            }
        }
    }
}


// Find maximum level of a shell.
void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const label shellI,
    labelList& maxLevel
) const
{
    const labelList& levels = levels_[shellI];

    if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.

        const scalarField& distances = distances_[shellI];

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(maxLevel, pointI)
        {
            forAllReverse(levels, levelI)
            {
                if (levels[levelI] > maxLevel[pointI])
                {
                    candidates[candidateI] = pt[pointI];
                    candidateMap[candidateI] = pointI;
                    candidateDistSqr[candidateI] = sqr(distances[levelI]);
                    candidateI++;
                    break;
                }
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        // Update maxLevel
        forAll(nearInfo, candidateI)
        {
            if (nearInfo[candidateI].hit())
            {
                // Check which level it actually is in.
                label minDistI = findLower
                (
                    distances,
                    mag(nearInfo[candidateI].hitPoint()-candidates[candidateI])
                );

                label pointI = candidateMap[candidateI];

                // pt is inbetween shell[minDistI] and shell[minDistI+1]
                maxLevel[pointI] = levels[minDistI+1];
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        label candidateI = 0;

        forAll(maxLevel, pointI)
        {
            if (levels[0] > maxLevel[pointI])
            {
                candidates[candidateI] = pt[pointI];
                candidateMap[candidateI] = pointI;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<searchableSurface::volumeType> volType;
        allGeometry_[shells_[shellI]].getVolumeType(candidates, volType);

        forAll(volType, i)
        {
            label pointI = candidateMap[i];

            if
            (
                (
                    modes_[shellI] == INSIDE
                 && volType[i] == searchableSurface::INSIDE
                )
             || (
                    modes_[shellI] == OUTSIDE
                 && volType[i] == searchableSurface::OUTSIDE
                )
            )
            {
                maxLevel[pointI] = levels[0];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shellSurfaces::shellSurfaces
(
    const searchableSurfaces& allGeometry,
    const PtrList<dictionary>& shellDicts
)
:
    allGeometry_(allGeometry)
{
    shells_.setSize(shellDicts.size());
    modes_.setSize(shellDicts.size());
    distances_.setSize(shellDicts.size());
    levels_.setSize(shellDicts.size());

    forAll(shellDicts, shellI)
    {
        const dictionary& dict = shellDicts[shellI];
        const word name = dict.lookup("name");
        const word type = dict.lookup("type");

        shells_[shellI] = allGeometry_.findSurfaceID(name);

        if (shells_[shellI] == -1)
        {
            FatalErrorIn
            (
                "shellSurfaces::shellSurfaces"
                "(const searchableSurfaces&, const PtrList<dictionary>&)"
            )   << "No surface called " << name << endl
                << "Valid surfaces are " << allGeometry_.names()
                << exit(FatalError);
        }

        modes_[shellI] = refineModeNames_.read(dict.lookup("refineMode"));

        // Read pairs of distance+level
        setAndCheckLevels(shellI, dict.lookup("levels"));
    }

    // Orient shell surfaces before any searching is done. Note that this
    // only needs to be done for inside or outside. Orienting surfaces
    // constructs lots of addressing which we want to avoid.
    orient();
}


Foam::shellSurfaces::shellSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& shellsDict
)
:
    allGeometry_(allGeometry)
{
    shells_.setSize(shellsDict.size());
    modes_.setSize(shellsDict.size());
    distances_.setSize(shellsDict.size());
    levels_.setSize(shellsDict.size());

    label shellI = 0;
    forAllConstIter(dictionary, shellsDict, iter)
    {
        shells_[shellI] = allGeometry_.findSurfaceID(iter().keyword());

        if (shells_[shellI] == -1)
        {
            FatalErrorIn
            (
                "shellSurfaces::shellSurfaces"
                "(const searchableSurfaces&, const dictionary>&"
            )   << "No surface called " << iter().keyword() << endl
                << "Valid surfaces are " << allGeometry_.names()
                << exit(FatalError);
        }
        const dictionary& dict = shellsDict.subDict(iter().keyword());

        modes_[shellI] = refineModeNames_.read(dict.lookup("mode"));

        // Read pairs of distance+level
        setAndCheckLevels(shellI, dict.lookup("levels"));

        shellI++;
    }

    // Orient shell surfaces before any searching is done. Note that this
    // only needs to be done for inside or outside. Orienting surfaces
    // constructs lots of addressing which we want to avoid.
    orient();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Highest shell level
Foam::label Foam::shellSurfaces::maxLevel() const
{
    label overallMax = 0;
    forAll(levels_, shellI)
    {
        overallMax = max(overallMax, max(levels_[shellI]));
    }
    return overallMax;
}


void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& maxLevel
) const
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;

    forAll(shells_, shellI)
    {
        findHigherLevel(pt, shellI, maxLevel);
    }
}


// ************************************************************************* //
