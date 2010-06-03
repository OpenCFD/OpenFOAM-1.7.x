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

Application
    surfaceFeatureExtract

Description
    Extracts and writes surface features to file

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "triSurface.H"
#include "argList.H"
#include "surfaceFeatures.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dumpBox(const treeBoundBox& bb, const fileName& fName)
{
    OFstream str(fName);

    Pout<< "Dumping bounding box " << bb << " as lines to obj file "
        << str.name() << endl;


    pointField boxPoints(bb.points());

    forAll(boxPoints, i)
    {
        meshTools::writeOBJ(str, boxPoints[i]);
    }

    forAll(treeBoundBox::edges, i)
    {
        const edge& e = treeBoundBox::edges[i];

        str<< "l " << e[0]+1 <<  ' ' << e[1]+1 << nl;
    }
}


// Deletes all edges inside/outside bounding box from set.
void deleteBox
(
    const triSurface& surf,
    const treeBoundBox& bb,
    const bool removeInside,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    forAll(edgeStat, edgeI)
    {
        const point eMid = surf.edges()[edgeI].centre(surf.localPoints());

        if
        (
            (removeInside && bb.contains(eMid))
         || (!removeInside && !bb.contains(eMid))
        )
        {
            edgeStat[edgeI] = surfaceFeatures::NONE;
        }
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

    argList::validArgs.clear();
    argList::validArgs.append("surface");
    argList::validArgs.append("output set");

    argList::validOptions.insert("includedAngle", "included angle [0..180]");
    argList::validOptions.insert("set", "input feature set");

    argList::validOptions.insert("minLen", "cumulative length of feature");
    argList::validOptions.insert("minElem", "number of edges in feature");
    argList::validOptions.insert("subsetBox", "((x0 y0 z0)(x1 y1 z1))");
    argList::validOptions.insert("deleteBox", "((x0 y0 z0)(x1 y1 z1))");
    argList args(argc, argv);

    Pout<< "Feature line extraction is only valid on closed manifold surfaces."
        << endl;


    fileName surfFileName(args.additionalArgs()[0]);
    fileName outFileName(args.additionalArgs()[1]);

    Pout<< "Surface            : " << surfFileName << nl
        << "Output feature set : " << outFileName << nl
        << endl;


    // Read
    // ~~~~

    triSurface surf(surfFileName);

    Pout<< "Statistics:" << endl;
    surf.writeStats(Pout);
    Pout<< endl;




    // Either construct features from surface&featureangle or read set.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    surfaceFeatures set(surf);

    if (args.optionFound("set"))
    {
        fileName setName(args.option("set"));

        Pout<< "Reading existing feature set from file " << setName << endl;

        set = surfaceFeatures(surf, setName);
    }
    else if (args.optionFound("includedAngle"))
    {
        scalar includedAngle = args.optionRead<scalar>("includedAngle");

        Pout<< "Constructing feature set from included angle " << includedAngle
            << endl;

        set = surfaceFeatures(surf, includedAngle);

        Pout<< endl << "Writing initial features" << endl;
        set.write("initial.fSet");
        set.writeObj("initial");
    }
    else
    {
        FatalErrorIn(args.executable())
            << "No initial feature set. Provide either one"
            << " of -set (to read existing set)" << nl
            << " or -includedAngle (to new set construct from angle)"
            << exit(FatalError);
    }


    Pout<< nl
        << "Initial feature set:" << nl
        << "    feature points : " << set.featurePoints().size() << nl
        << "    feature edges  : " << set.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << set.nRegionEdges() << nl
        << "        external edges : " << set.nExternalEdges() << nl
        << "        internal edges : " << set.nInternalEdges() << nl
        << endl;




    // Trim set
    // ~~~~~~~~

    scalar minLen = -GREAT;
    if (args.optionReadIfPresent("minLen", minLen))
    {
        Pout<< "Removing features of length < " << minLen << endl;
    }

    label minElem = 0;
    if (args.optionReadIfPresent("minElem", minElem))
    {
        Pout<< "Removing features with number of edges < " << minElem << endl;
    }

    // Trim away small groups of features
    if (minLen > 0 || minLen > 0)
    {
        set.trimFeatures(minLen, minElem);
        Pout<< endl << "Removed small features" << endl;
    }



    // Subset
    // ~~~~~~

    // Convert to marked edges, points
    List<surfaceFeatures::edgeStatus> edgeStat(set.toStatus());

    if (args.optionFound("subsetBox"))
    {
        treeBoundBox bb(args.optionLookup("subsetBox")());

        Pout<< "Removing all edges outside bb " << bb << endl;
        dumpBox(bb, "subsetBox.obj");

        deleteBox
        (
            surf,
            bb,
            false,
            edgeStat
        );
    }
    else if (args.optionFound("deleteBox"))
    {
        treeBoundBox bb(args.optionLookup("deleteBox")());

        Pout<< "Removing all edges inside bb " << bb << endl;
        dumpBox(bb, "deleteBox.obj");

        deleteBox
        (
            surf,
            bb,
            true,
            edgeStat
        );
    }

    surfaceFeatures newSet(surf);
    newSet.setFromStatus(edgeStat);

    Pout<< endl << "Writing trimmed features to " << outFileName << endl;
    newSet.write(outFileName);

    Pout<< endl << "Writing edge objs." << endl;
    newSet.writeObj("final");


    Pout<< nl
        << "Final feature set:" << nl
        << "    feature points : " << newSet.featurePoints().size() << nl
        << "    feature edges  : " << newSet.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << newSet.nRegionEdges() << nl
        << "        external edges : " << newSet.nExternalEdges() << nl
        << "        internal edges : " << newSet.nInternalEdges() << nl
        << endl;

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
