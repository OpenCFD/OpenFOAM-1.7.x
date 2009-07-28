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

Application
    foamToStarMesh

Description
    Reads an OpenFOAM mesh and writes a pro-STAR (v4) bnd/cel/vrt format.

Usage
    - foamToStarMesh [OPTION] \n
    Reads an OpenFOAM mesh and writes a pro-STAR (v4) bnd/cel/vrt format.

    @param -noBnd \n
    Suppress writing the @c .bnd file

    @param -scale \<factor\>\n
    Specify an alternative geometry scaling factor.
    The default is @b 1000 (scale @em [m] to @em [mm]).

    @param -surface \n
    Extract the surface of the volume mesh only.
    This can be useful, for example, for surface morphing in an external
    package.

    @param -tri \n
    Extract a triangulated surface.
    The @b -surface options is implicitly selected.


Note
    The cellTable information available in the files
    @c constant/cellTable and @c constant/polyMesh/cellTableId
    will be used if available. Otherwise the cellZones are used when
    creating the cellTable information.

See Also
    Foam::cellTable, Foam::meshWriter and Foam::meshWriters::STARCD

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "STARCDMeshWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

    argList::validOptions.insert("scale", "scale");
    argList::validOptions.insert("noBnd", "");
    argList::validOptions.insert("tri", "");
    argList::validOptions.insert("surface", "");

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    bool surfaceOnly = false;
    if (args.optionFound("surface") || args.optionFound("tri"))
    {
        surfaceOnly = true;
    }

    fileName exportName = meshWriter::defaultMeshName;
    if (surfaceOnly)
    {
        exportName = meshWriter::defaultSurfaceName;
    }

    if (args.optionFound("case"))
    {
        exportName += '-' + args.globalCaseName();
    }

    // default: rescale from [m] to [mm]
    scalar scaleFactor = 1000;
    if (args.optionReadIfPresent("scale", scaleFactor))
    {
        if (scaleFactor <= 0)
        {
            scaleFactor = 1;
        }
    }

#   include "createPolyMesh.H"


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

#       include "getTimeIndex.H"

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != polyMesh::UNCHANGED)
        {
            meshWriters::STARCD writer(mesh, scaleFactor);

            if (args.optionFound("noBnd"))
            {
                writer.noBoundary();
            }

            fileName meshName(exportName);
            if (state != polyMesh::UNCHANGED)
            {
                meshName += '_' + runTime.timeName();
            }

            if (surfaceOnly)
            {
                if (args.optionFound("tri"))
                {
                    writer.writeSurface(meshName, true);
                }
                else
                {
                    writer.writeSurface(meshName);
                }
            }
            else
            {
                writer.write(meshName);
            }
        }

        Info<< nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
