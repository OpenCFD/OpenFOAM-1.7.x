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

#include "VTKsurfaceFormat.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::writeHeaderPolygons
(
    Ostream& os,
    const UList<Face>& faceLst
)
{
    label nNodes = 0;

    forAll(faceLst, faceI)
    {
        nNodes += faceLst[faceI].size();
    }

    os  << nl
        << "POLYGONS " << faceLst.size() << ' '
        << faceLst.size() + nNodes << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::VTKsurfaceFormat<Face>::VTKsurfaceFormat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>&  faceLst = surf.faces();
    const List<label>& faceMap = surf.faceMap();

    const List<surfZone>& zones =
    (
        surf.surfZones().size() > 1
      ? surf.surfZones()
      : oneZone(faceLst)
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorIn
        (
            "fileFormats::VTKsurfaceFormat::write"
            "(const fileName&, const MeshedSurfaceProxy<Face>&)"
        )
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    writeHeader(os, pointLst);
    writeHeaderPolygons(os, faceLst);

    label faceIndex = 0;
    forAll(zones, zoneI)
    {
        const surfZone& zone = zones[zoneI];

        if (useFaceMap)
        {
            forAll(zone, localFaceI)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                os << f.size();
                forAll(f, fp)
                {
                    os << ' ' << f[fp];
                }
                os << ' ' << nl;
            }
        }
        else
        {
            forAll(zone, localFaceI)
            {
                const Face& f = faceLst[faceIndex++];

                os << f.size();
                forAll(f, fp)
                {
                    os << ' ' << f[fp];
                }
                os << ' ' << nl;
            }
        }
    }

    writeTail(os, zones);
}


template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorIn
        (
            "fileFormats::VTKsurfaceFormat::write"
            "(const fileName&, const UnsortedMeshedSurface<Face>&)"
        )
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    const List<Face>& faceLst = surf.faces();

    writeHeader(os, surf.points());
    writeHeaderPolygons(os, faceLst);

    forAll(faceLst, faceI)
    {
        const Face& f = faceLst[faceI];

        os << f.size();
        forAll(f, fp)
        {
            os << ' ' << f[fp];
        }
        os << ' ' << nl;
    }

    writeTail(os, surf.zoneIds());
}


// ************************************************************************* //
