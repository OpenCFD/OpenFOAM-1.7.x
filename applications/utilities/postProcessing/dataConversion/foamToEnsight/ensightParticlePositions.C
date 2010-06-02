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

#include "ensightParticlePositions.H"
#include "fvMesh.H"
#include "passiveParticle.H"
#include "Cloud.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "itoa.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void ensightParticlePositions
(
    const Foam::fvMesh& mesh,
    const Foam::fileName& postProcPath,
    const Foam::word& timeFile,
    const Foam::word& cloudName,
    const bool dataExists
)
{
    if (dataExists)
    {
        Info<< "Converting cloud " << cloudName << " positions" <<  endl;
    }
    else
    {
        Info<< "Creating empty cloud " << cloudName << " positions" << endl;
    }

    const Time& runTime = mesh.time();

    fileName ensightFileName(timeFile + "." + cloudName);
    OFstream ensightFile
    (
        postProcPath/ensightFileName,
        runTime.writeFormat(),
        runTime.writeVersion(),
        runTime.writeCompression()
    );

    // Output header
    ensightFile
        << cloudName.c_str() << nl
        << "particle coordinates" << nl;

    if (dataExists)
    {
        Cloud<passiveParticle> parcels(mesh, cloudName, false);

        // Set Format
        ensightFile.setf(ios_base::scientific, ios_base::floatfield);
        ensightFile.precision(5);

        ensightFile<< setw(8) << parcels.size() << nl;

        label nParcels = 0;

        // Output positions
        forAllConstIter(Cloud<passiveParticle>, parcels, elmnt)
        {
            const vector& p = elmnt().position();

            ensightFile
                << setw(8) << ++nParcels
                << setw(12) << p.x() << setw(12) << p.y() << setw(12) << p.z()
                << nl;
        }
    }
    else
    {
        label nParcels = 0;
        ensightFile<< setw(8) << nParcels << nl;
    }
}


// ************************************************************************* //
