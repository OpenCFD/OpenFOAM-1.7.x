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

#include "STARCDsurfaceFormatCore.H"
#include "clock.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileFormats::STARCDsurfaceFormatCore::readHeader
(
    IFstream& is,
    const word& signature
)
{
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::STARCDsurfaceFormatCore::readHeader(...)"
        )
            << "cannot read " << signature  << "  " << is.name()
            << abort(FatalError);
    }

    word header;
    label majorVersion;

    string line;

    is.getLine(line);
    IStringStream(line)() >> header;

    is.getLine(line);
    IStringStream(line)() >> majorVersion;

    // add other checks ...
    if (header != signature)
    {
        Info<< "header mismatch " << signature << "  " << is.name()
            << endl;
    }

    return true;
}


void Foam::fileFormats::STARCDsurfaceFormatCore::writeHeader
(
    Ostream& os,
    const char* filetype
)
{
    os  << "PROSTAR_" << filetype << nl
        << 4000
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << endl;
}


bool Foam::fileFormats::STARCDsurfaceFormatCore::readPoints
(
    IFstream& is,
    pointField& points,
    labelList& ids
)
{
    //
    // read .vrt file
    // ~~~~~~~~~~~~~~

    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::STARCDsurfaceFormatCore::readPoints(...)"
        )
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    readHeader(is, "PROSTAR_VERTEX");

    DynamicList<point> dynPoints;
    // STAR-CD index of points
    DynamicList<label> dynPointId;

    label lineLabel;
    while ((is >> lineLabel).good())
    {
        scalar x, y, z;

        is >> x >> y >> z;

        dynPoints.append(point(x, y, z));
        dynPointId.append(lineLabel);
    }

    points.transfer(dynPoints);
    ids.transfer(dynPointId);

    return true;
}



void Foam::fileFormats::STARCDsurfaceFormatCore::writePoints
(
    Ostream& os,
    const pointField& pointLst
)
{
    writeHeader(os, "VERTEX");

    // Set the precision of the points data to 10
    os.precision(10);

    // force decimal point for Fortran input
    os.setf(std::ios::showpoint);

    forAll(pointLst, ptI)
    {
        os
            << ptI + 1 << " "
            << pointLst[ptI].x() << " "
            << pointLst[ptI].y() << " "
            << pointLst[ptI].z() << nl;
    }
    os.flush();
}


void Foam::fileFormats::STARCDsurfaceFormatCore::writeCase
(
    Ostream& os,
    const pointField& pointLst,
    const label nFaces,
    const UList<surfZone>& zoneLst
)
{
    word caseName = os.name().lessExt().name();

    os  << "! STAR-CD file written " << clock::dateTime().c_str() << nl
        << "! " << pointLst.size() << " points, " << nFaces << " faces" << nl
        << "! case " << caseName << nl
        << "! ------------------------------" << nl;

    forAll(zoneLst, zoneI)
    {
        os  << "ctable " << zoneI + 1 << " shell" << nl
            << "ctname " << zoneI + 1 << " "
            << zoneLst[zoneI].name() << nl;
    }

    os  << "! ------------------------------" << nl
        << "*set icvo mxv - 1" << nl
        << "vread " << caseName << ".vrt icvo,,,coded" << nl
        << "cread " << caseName << ".cel icvo,,,add,coded" << nl
        << "*set icvo" << nl
        << "! end" << nl;

    os.flush();
}


// ************************************************************************* //

