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
    surfaceFeatureConvert

Description
    Extracts and writes surface features to file

\*---------------------------------------------------------------------------*/

#include "featureEdgeMesh.H"
#include "argList.H"
#include "Time.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "OFstream.H"
#include "Map.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void readNASEdges
(
    const fileName& inFileName,
    pointField& allPoints,
    edgeList& allEdges
)
{
    IFstream is(inFileName);

    if (!is.good())
    {
        FatalErrorIn("readNASEdges")
            << "Cannot read file " << inFileName
            << exit(FatalError);
    }

    // coordinates of point
    DynamicList<point> points;
    // Nastran index of point
    DynamicList<label> pointIndices;

    // beams
    DynamicList<edge> edges;
    DynamicList<label> edgeIndices;


    while (is.good())
    {
        string line;
        is.getLine(line);

        if (line.empty() || line[0] == '$')
        {
            // Skip empty and comment
            continue;
        }

        // Check if character 72 is continuation
        if (line.size() > 72 && line[72] == '+')
        {
            line = line.substr(0, 72);

            while (true)
            {
                string buf;
                is.getLine(buf);

                if (buf.size() > 72 && buf[72] == '+')
                {
                    line += buf.substr(8, 64);
                }
                else
                {
                    line += buf.substr(8, buf.size()-8);
                    break;
                }
            }
        }

        // Read first word
        IStringStream lineStream(line);
        word cmd;
        lineStream >> cmd;

        if (cmd == "GRID")
        {
            label index;
            lineStream >> index;
            pointIndices.append(index);

            scalar x = readScalar(IStringStream(line.substr(24, 8))());
            scalar y = readScalar(IStringStream(line.substr(32, 8))());
            scalar z = readScalar(IStringStream(line.substr(40, 8))());
            points.append(point(x, y, z));
        }
        else if (cmd == "CBEAM")
        {
            // Read shell type since gives patchnames.
            label index, group, v0, v1;
            lineStream >> index >> group >> v0 >> v1;

            edgeIndices.append(index);
            edges.append(edge(v0, v1));
        }
    }

    points.shrink();
    pointIndices.shrink();
    edges.shrink();
    edgeIndices.shrink();

    Pout<< "Read from " << inFileName
        << " edges:" << edges.size() << " points:" << points.size()
        << endl;

    {
        // Build inverse mapping (index to point)
        Map<label> indexToPoint(2*pointIndices.size());
        forAll(pointIndices, i)
        {
            indexToPoint.insert(pointIndices[i], i);
        }

        // Relabel edges
        forAll(edges, i)
        {
            edge& e = edges[i];
            e[0] = indexToPoint[e[0]];
            e[1] = indexToPoint[e[1]];
        }
    }

    allPoints.transfer(points);
    allEdges.transfer(edges);
}



void write
(
    const Time& runTime,
    const fileName& inFileName,
    const fileName& outFileName,
    const edgeMesh& eMesh
)
{
    if (outFileName.ext() == "eMesh")
    {
        featureEdgeMesh fem
        (
            IOobject
            (
                outFileName,        // name
                runTime.constant(), // instance
                runTime,            // registry
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            eMesh.points(),
            eMesh.edges()
        );

        Pout<< "Writing feature edge mesh to " << fem.objectPath()
            << endl;

        fem.write();
    }
    else if (outFileName.ext() == "vtk")
    {
        OFstream str(outFileName);

        str << "# vtk DataFile Version 2.0" << nl
            << "featureEdgeMesh " << inFileName << nl
            << "ASCII" << nl
            << "DATASET POLYDATA" << nl;

        str << "POINTS " << eMesh.points().size() << " float" << nl;
        forAll(eMesh.points(), pointI)
        {
            const point& pt = eMesh.points()[pointI];

            str << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
        }

        str << "LINES " << eMesh.edges().size() << ' '
            << 3*eMesh.edges().size() << nl;
        forAll(eMesh.edges(), edgeI)
        {
            const edge& e = eMesh.edges()[edgeI];

            str << "2 " << e[0] << ' ' << e[1] << nl;
        }
    }
    else
    {
        FatalErrorIn("write")
            << "Supported output formats: .eMesh, .vtk"
            << exit(FatalError);
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("input file");
    argList::validArgs.append("output file");
#   include "setRootCase.H"
#   include "createTime.H"

    const fileName inFileName(args.additionalArgs()[0]);
    const word outFileName(args.additionalArgs()[1]);

    Pout<< "Input features file  : " << inFileName << nl
        << "Output features file : " << outFileName << nl
        << endl;


    // Read
    // ~~~~

    if (inFileName.ext() == "nas")
    {
        pointField points;
        edgeList edges;
        readNASEdges(inFileName, points, edges);

        edgeMesh eMesh(points, edges);

        write(runTime, inFileName, outFileName, eMesh);
    }
    else if (inFileName.ext() == "eMesh")
    {
        featureEdgeMesh fem
        (
            IOobject
            (
                inFileName,         // name
                runTime.constant(), // instance
                runTime,            // registry
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            )
        );

        Pout<< "Read from " << inFileName
            << " edges:" << fem.edges().size()
            << " points:" << fem.points().size()
            << endl;

        write(runTime, inFileName, outFileName, fem);
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Can only handle NASTRAN data formats (.nas extension)."
            << exit(FatalError);
    }

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
